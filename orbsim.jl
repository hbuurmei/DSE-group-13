using Plots:display
using Base:print_array
using Pkg
# Pkg.add("LoopVectorization")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")
# Pkg.add("PyPlot")

using CSV
using DataFrames
using Statistics
using LinearAlgebra
using LoopVectorization
using Plots

# Fundamental constants
const R_e = 6378.137e3  # [m]
const g_0 = 9.80665  # [m/s2]
const J_2 = 0.00108263  # [-]
const mu = 3.986004418e14  # [m3/s2]

# User defined constants
const h_collision = 789e3  # [m]
const debris_n = 1000  # number of fragments, change this number for simulation speed
const a_collision = R_e + h_collision
const t0 = 72 * 100 * 60  # 5 days
const dt = 5
const distance_sc = 30e3  # [m]
const target_fraction = 0.5
const max_dv = 1 # Maximum dV used in gaussian perturbation equations
const FoV = 54.63 * pi / 180  # [rad]
const range = 250e3 # [m]
const incidence_angle = 20 * pi / 180 # [rad]
const ablation_time = 50 # [s]
const scan_time = 10 # [s]
const min_vis_time = scan_time + ablation_time # [s]
const cooldown_time = min_vis_time + 0 # seconds, should be an integer multiple of dt
const view_angles = (45, 45) # Viewing angles in azimuth and altitude

# Import data
df = CSV.read("./iridium_cosmos_result.csv", DataFrame; header=1)

# Select data that is necessary and convert to matrix
df = filter(row -> row.Name .== "Kosmos 2251-Collision-Fragment", df)
df = filter(row -> row.d_eq .< 0.1, df)
df = filter(row -> 0 .< row.e .< 1, df)
df = filter(row -> (row.a * (1 - row.e) .> (R_e + 200e3)) && (row.a * (1 + row.e) .> (R_e + 200e3)), df) # Filter out all that already have a low enough perigee
debris_kepler = Matrix(select(df, ["a", "e", "i", "long_asc", "arg_peri", "mean_anom", "ID"])) # ID is used as an additional column to store true anomaly
debris_dims = Matrix(select(df, ["M"]))

# Cut data set down to set number of fragments
tot_debris_n = min(debris_n, length(debris_kepler[:,1]))
println("Number of debris objects: ", tot_debris_n)
debris_kepler = debris_kepler[1:tot_debris_n,:]
debris_semimajor_original = debris_kepler[:,1]
debris_cartesian = Matrix{Float64}(undef, tot_debris_n, 3)
debris_cartesian_vel = Matrix{Float64}(undef, tot_debris_n, 3)
debris_removed = zeros(Bool, tot_debris_n, 2)
debris_vis_times_pass = zeros(tot_debris_n)
debris_vis_prev = zeros(Bool, tot_debris_n)

# Spacecraft constants
const a_sc = R_e + h_collision + distance_sc
const e_sc = 0
const laser_pointing_angle = acos(a_collision / a_sc)  # [rad]
const M_0_sc = 0  # Assume SC won't be exactly in phase with debris
const i_sc = 74 * pi / 180 # [rad], Kosmos inclination


@inline function calc_true_anomaly(a, e, M)
    # Initial guess
    E = 0

    # Apply newton method 5x (reaches max precision after 5 iterations)
    for i = 1:5
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E))
    end

    # Final equation for true anomaly
    return 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
end

# TODO Perhaps pass the array into function and just fill directly
function kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN, position)
    # Convert a position in the Keplerian system to a cartesian system
    p = a * (1 - e * e)
    r = p / (1 + e * cos(true_anomaly)) # radius

    # Compute the Cartesian position vector
    X = r * (cos(RAAN) * cos(w + true_anomaly) - sin(RAAN) * sin(w + true_anomaly) * cos(i))
    Y = r * (sin(RAAN) * cos(w + true_anomaly) + cos(RAAN) * sin(w + true_anomaly) * cos(i))
    Z = r * (sin(i) * sin(w + true_anomaly))
    
    position[1] = X
    position[2] = Y
    position[3] = Z
end

function calc_vel(a, e, w, true_anomaly, i, RAAN, position)
    # Get the velocity in cartesian coordinates
    p = a * (1 - e * e)
    r = p / (1 + e * cos(true_anomaly)) # radius
    h = sqrt(mu * p)

    V_X = (position[1] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (cos(RAAN) * sin(w + true_anomaly) + sin(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Y = (position[2] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (sin(RAAN) * sin(w + true_anomaly) - cos(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Z = (position[3] * h * e / (r * p)) * sin(true_anomaly) + (h / r) * (cos(w + true_anomaly) * sin(i))

    return [V_X, V_Y, V_Z]
end

function thrust_alter_orbit(debris_kepler, debris_cartesian, debris_cartesian_vel, debris_dims, thrust_dir, thrust_energy, i)
    # Establish RTO (Radial, Transverse, Out-of-plane) axes (unit vectors)
    @inbounds R = normalize(debris_cartesian[i,:])
    @inbounds O = normalize(cross(R, debris_cartesian_vel[i,:]))
    T = cross(O, R)
    
    # Compute product of a and thrust_dt
    # Based on kinetic energy and v2 = a * dt + v1
    @inbounds v1 = norm(debris_cartesian_vel[i,:])
    @inbounds tot_dv = sqrt(v1 * v1 + 2 * thrust_energy / debris_dims[i,1]) - v1

    remaining_dv = tot_dv
    while remaining_dv > 0
        dv = (remaining_dv / max_dv) < 1 ? mod(remaining_dv, max_dv) : max_dv
        dir_dv = normalize(thrust_dir) .* dv
        dir_dv_rto = zeros(3)
        @inbounds dir_dv_rto[1] = dot(dir_dv, R)
        @inbounds dir_dv_rto[2] = dot(dir_dv, T)
        @inbounds dir_dv_rto[3] = dot(dir_dv, O)

        @inbounds sqramu = sqrt(debris_kepler[i, 1] / mu)
        @inbounds sub1e2 = 1 - debris_kepler[i, 2] * debris_kepler[i, 2]
        sqr1e2 = sqrt(sub1e2)
        @inbounds sinf = sin(debris_kepler[i, 7])
        @inbounds cosf = cos(debris_kepler[i, 7])
        @inbounds ecosf1 = debris_kepler[i, 2] * cosf + 1
        @inbounds n = sqrt(mu / debris_kepler[i, 1]^3)

        # Gaussian perturbation formulae
        @inbounds debris_kepler[i, 1] += sqramu * 2 * debris_kepler[i, 1] / sqr1e2 * (debris_kepler[i, 2] * sinf * dir_dv_rto[1] + ecosf1 * dir_dv_rto[2])
        @inbounds debris_kepler[i, 2] += sqramu * sqr1e2 * (sinf * dir_dv_rto[1] + (debris_kepler[i, 2] + 2 * cosf + debris_kepler[i, 2] * cosf * cosf) / ecosf1 * dir_dv_rto[2])
        @inbounds debris_kepler[i, 3] += sqramu * sqr1e2 / ecosf1 * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) * dir_dv_rto[3]
        dRAAN = sqramu * sqr1e2 / ecosf1 * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) / sin(debris_kepler[i, 3]) * dir_dv_rto[3]
        @inbounds debris_kepler[i, 4] += dRAAN
        @inbounds debris_kepler[i, 5] += sqramu * sqr1e2 / debris_kepler[i, 2] * (- cosf * dir_dv_rto[1] + (ecosf1 + 1) / ecosf1 * sinf * dir_dv_rto[2]) - cos(debris_kepler[i, 3]) * dRAAN
        @inbounds debris_kepler[i, 6] += n + sub1e2 / (n * debris_kepler[i, 1] * debris_kepler[i, 2]) * ((cosf - 2 * debris_kepler[i, 2] / ecosf1) * dir_dv_rto[1] - (ecosf1 + 1) / ecosf1 * sinf * dir_dv_rto[2])

        remaining_dv -= max_dv
    end

    # println("ΔV imparted: ", (sqrt(v1 * v1 + 2 * thrust_energy / debris_dims[i,1]) - v1))
end

function J_2_RAAN(a, e, i)
    n = sqrt(mu / a^3)
    RAAN_dot = -1.5 * n * R_e * R_e * J_2 * cos(i) / (a * a) / (1 - e * e)^2
    return RAAN_dot
end

function J_2_w(a, e, i)
    n = sqrt(mu / a^3)
    w_dot = 0.75 * n * R_e * R_e * J_2 * (4 - 5 * (sin(i))^2) / (a * a) / (1 - e * e)^2
    return w_dot
end

function run_sim(;plotResults=true)
    debris_counter = 0
    increased_a_counter = 0
    t = t0
    t_last_pulse = -Inf64
    w_sc = J_2_w(a_sc, e_sc, i_sc) * t0
    # RAAN_sc = mean(debris_kepler[:, 4]) + J_2_RAAN(a_sc, e_sc, i_sc) * t0
    ts = Vector{Float64}(undef, 0)
    percentages = Vector{Float64}(undef, 0)
    position_sc = zeros(3)
    vel_sc = zeros(3)
    camera_axis_dot = zeros(tot_debris_n)

    temp_debris_kepler = debris_kepler[1,:]

    sizehint!(ts, 100000);
    sizehint!(percentages, 100000);

    # J_2 effect sc
    RAAN_drift_sc = J_2_RAAN(a_sc, e_sc, i_sc) * dt
    w_drift_sc = J_2_w(a_sc, e_sc, i_sc) * dt

    # J_2 effect debris
    RAAN_drift = Vector{Float64}(undef, tot_debris_n)
    w_drift = Vector{Float64}(undef, tot_debris_n)

    # Precompute RAAN and w drifts due to J2 effect
    for i in eachindex(RAAN_drift)
        @inbounds RAANd = J_2_RAAN(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3])
        @inbounds wd = J_2_w(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3])
        @inbounds RAAN_drift[i] = RAANd * dt
        @inbounds w_drift[i] = wd * dt
        @inbounds debris_kepler[i, 4] += RAANd * t0
        @inbounds debris_kepler[i, 5] += wd * t0
    end

    # Assume SC can be inserted at similar RAAN angle as space debris
    RAAN_sc = mean(debris_kepler[:, 4])

    while (debris_counter / tot_debris_n < target_fraction)
        push!(ts, t - t0)

        # Update RAAN and w due to J_2 (sc)
        RAAN_sc += RAAN_drift_sc
        w_sc += w_drift_sc

        # Compute spacecraft position and velocity
        n_sc = sqrt(mu / a_sc^3)
        true_anomaly_sc = calc_true_anomaly(a_sc, e_sc, n_sc * t + M_0_sc)
        kepler_to_cartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)
        vel_sc = calc_vel(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)

        # Update space debris position
        @turbo for i = 1:tot_debris_n
            # left here for readability
            # a = debris_kepler[i, 1], semi-major axis
            # e = debris_kepler[i, 2], eccentricity
            # inc = debris_kepler[i, 3], inclination
            # RAAN = debris_kepler[i, 4], right ascension of ascending node
            # w = debris_kepler[i, 5], argument of pericenter
            # M = debris_kepler[i, 6], mean anomaly
            # f = debris_kepler[i, 7], true anomaly

            # Update RAAN and w due to J_2 (debris)
            @inbounds debris_kepler[i, 4] += RAAN_drift[i]
            @inbounds debris_kepler[i, 5] += w_drift[i]

            # Update mean anomaly
            @inbounds n = sqrt(mu / debris_kepler[i, 1]^3)
            @inbounds debris_kepler[i, 6] = mod(debris_kepler[i, 6] + n * dt, 2 * pi)

            # Update true anomaly
            @inbounds debris_kepler[i, 7] = calc_true_anomaly(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 6])

            @inbounds p = debris_kepler[i, 1] * (1 - debris_kepler[i, 2] * debris_kepler[i, 2])
            @inbounds r = p / (1 + debris_kepler[i, 2] * cos(debris_kepler[i, 7])) # radius

            # Compute the Cartesian position vector
            @inbounds debris_cartesian[i, 1] = r * (cos(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) - sin(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 2] = r * (sin(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) + cos(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 3] = r * (sin(debris_kepler[i, 3]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]))
        end

        # This is separate from the above loop because @tturbo uses vector intrinsics, which are not available for more complex functions
        for i = 1:tot_debris_n
            
            if t - t_last_pulse < cooldown_time
                break # If laser is not ready, skip rest of the loop
            elseif debris_removed[i,1]
                continue # If debris object is already marked as removed, skip it
            end

            @inbounds rel_pos = position_sc - debris_cartesian[i,:] # Vector from debris to spacecraft
            @inbounds abs_distance = norm(rel_pos)
            if abs_distance < range
                # Update spacecraft velocity
                @inbounds debris_cartesian_vel[i,:] = calc_vel(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 5], debris_kepler[i, 7], debris_kepler[i, 3], debris_kepler[i, 4], debris_cartesian[i,:])
                
                # Check angle between debris tranjectory and spacecraft relative to debris
                # println(dot(debris_velocity, rel_pos) / (norm(debris_velocity) * norm(rel_pos)))
                @inbounds vel_rel_pos_angle = acos(sum(debris_cartesian_vel[i,:] .* rel_pos) / (norm(debris_cartesian_vel[i,:]) * norm(rel_pos)))
                incidence_condition = vel_rel_pos_angle < incidence_angle
                # @inbounds angvel_condition = (norm(debris_cartesian_vel[i,:]) / abs_distance * sin(vel_rel_pos_angle)) < 2 * pi / 180
                rotation_vector = cross(position_sc, -vel_sc) / norm(cross(position_sc, -vel_sc))
                pointing_vector = -vel_sc .* cos(laser_pointing_angle) + cross(rotation_vector, -vel_sc) .* sin(laser_pointing_angle) + rotation_vector .* (dot(rotation_vector, -vel_sc) * (1 - cos(laser_pointing_angle)))
                FoV_condition = acos(dot(pointing_vector, -rel_pos) / (norm(pointing_vector) * norm(-rel_pos))) < FoV / 2
                if incidence_condition && FoV_condition
                    # Inside sphere and cone
                    # println("Inside cone")
                    

                    if !debris_vis_prev[i]
                        temp_position_sc = position_sc
                        temp_vel_sc = vel_sc
                        temp_w_sc = w_sc
                        temp_RAAN_sc = RAAN_sc
                        
                        @inbounds temp_debris_kepler = debris_kepler[i,:]
                        @inbounds temp_debris_cartesian = debris_cartesian[i,:]
                        @inbounds temp_debris_cartesian_vel = debris_cartesian_vel[i,:]
                        temp_rel_pos = position_sc
                        predicted_vis_time = 0

                        temp_indicence_condition = true
                        temp_FoV_condition = true
                        temp_range_condition = true
                        while temp_range_condition && temp_indicence_condition# && temp_FoV_condition
                            predicted_vis_time += dt

                            # Propagate spacecraft object forward once
                            temp_w_sc += w_drift_sc
                            temp_RAAN_sc += RAAN_drift_sc
                            temp_n_sc = sqrt(mu / a_sc^3)
                            temp_true_anomaly_sc = calc_true_anomaly(a_sc, e_sc, n_sc * t + M_0_sc)
                            kepler_to_cartesian(a_sc, e_sc, temp_w_sc, temp_true_anomaly_sc, i_sc, temp_RAAN_sc, temp_position_sc)
                            temp_vel_sc = calc_vel(a_sc, e_sc, temp_w_sc, temp_true_anomaly_sc, i_sc, temp_RAAN_sc, temp_position_sc)

                            # Propagate debris object forward once
                            # Update RAAN and w due to J_2 (debris)
                            @inbounds temp_debris_kepler[4] += RAAN_drift[i]
                            @inbounds temp_debris_kepler[5] += w_drift[i]
                            # Update mean anomaly
                            @inbounds n = sqrt(mu / temp_debris_kepler[1]^3)
                            @inbounds temp_debris_kepler[6] = mod(temp_debris_kepler[6] + n * dt, 2 * pi)
                            # Update true anomaly
                            @inbounds temp_debris_kepler[7] = calc_true_anomaly(temp_debris_kepler[1], temp_debris_kepler[2], temp_debris_kepler[6])
                            # Update cartesian position
                            @inbounds kepler_to_cartesian(temp_debris_kepler[1], temp_debris_kepler[2], temp_debris_kepler[5], temp_debris_kepler[7], temp_debris_kepler[3], temp_debris_kepler[4], temp_debris_cartesian)
                            # Update cartesian velocity
                            @inbounds temp_debris_cartesian_vel = calc_vel(temp_debris_kepler[1], temp_debris_kepler[2], temp_debris_kepler[5], temp_debris_kepler[7], temp_debris_kepler[3], temp_debris_kepler[4], temp_debris_cartesian)

                            # Determine new relative position and distance
                            temp_rel_pos = position_sc - temp_debris_cartesian # Vector from debris to spacecraft
                            temp_abs_distance = norm(temp_rel_pos)

                            # Determine if range condition for new position is met
                            temp_range_condition = temp_abs_distance < range

                            # Determine if incidence condition for new position is met
                            temp_rel_pos = position_sc - temp_debris_cartesian
                            temp_vel_rel_pos_angle = acos(sum(temp_debris_cartesian_vel .* temp_rel_pos) / (norm(temp_debris_cartesian_vel) * norm(temp_rel_pos)))
                            temp_incidence_condition = temp_vel_rel_pos_angle < incidence_angle

                            # Determine if FoV condition for new position is met
                            temp_rotation_vector = cross(temp_position_sc, - temp_vel_sc) / norm(cross(temp_position_sc, - temp_vel_sc))
                            temp_pointing_vector = - temp_vel_sc .* cos(laser_pointing_angle) + cross(temp_rotation_vector, - temp_vel_sc) .* sin(laser_pointing_angle) + temp_rotation_vector .* (dot(temp_rotation_vector,  - temp_vel_sc) * (1 - cos(laser_pointing_angle)))
                            temp_FoV_condition = acos(dot(temp_pointing_vector, - temp_rel_pos) / (norm(temp_pointing_vector) * norm(temp_rel_pos))) < FoV / 2
                        end
                        debris_vis_times_pass[i] = predicted_vis_time
                        println("Fragment detected, expected: ", debris_vis_times_pass[i], " s in view.")
                    end

                    if debris_vis_times_pass[i] >= min_vis_time
                        @inbounds debris_removed[i,2] = true

                        @inbounds thrust_dir = - normalize(debris_cartesian_vel[i,:]) # Thrust opposite of debris velocity
                        energy_per_pulse = 5000 # J

                        @inbounds curr_true_anom = debris_kepler[i, 7] * 180 / pi
                        @inbounds curr_alt = (debris_kepler[i, 1] * (1 - debris_kepler[i, 2] * debris_kepler[i, 2]) / (1 + debris_kepler[i, 2] * cos(debris_kepler[i, 7])) - R_e)
                        @inbounds prev_perigee_alt = (debris_kepler[i, 1] * (1 - debris_kepler[i, 2]) - R_e)
                        @inbounds prev_apogee_alt = (debris_kepler[i, 1] * (1 + debris_kepler[i, 2]) - R_e)
                        thrust_alter_orbit(debris_kepler, debris_cartesian, debris_cartesian_vel, debris_dims, thrust_dir, energy_per_pulse, i)
                        @inbounds new_perigee_alt = (debris_kepler[i, 1] * (1 - debris_kepler[i, 2]) - R_e)
                        @inbounds new_apogee_alt = (debris_kepler[i, 1] * (1 + debris_kepler[i, 2]) - R_e)


                        # Update drifts
                        @inbounds RAAN_drift[i] = J_2_RAAN(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * dt
                        @inbounds w_drift[i] = J_2_w(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * dt

                        # println("Current True Anomaly: ", round(curr_true_anom, digits=0),"[deg], Current alt: ", round(curr_alt/1000, digits=2), "[km]")
                        # println("Old perigree alt: ", round(prev_perigee_alt/1000, digits=2), "[km], New perigee alt: ", round(new_perigee_alt/1000, digits=2), "[km]")
                        # println("Old apogree alt: ", round(prev_apogee_alt/1000, digits=2), "[km], New apogee alt: ", round(new_apogee_alt/1000, digits=2), "[km]")
                        @inbounds debris_removed[i,1] = (new_perigee_alt < 200e3) || (new_apogee_alt < 200e3) # Mark object as removed if perigee is now below 200 km
                        @inbounds debris_counter += debris_removed[i,1]
                        @inbounds increased_a_counter += (debris_semimajor_original[i] > a_collision)
                        println("New Perigee: ", new_perigee_alt, ", New Apogee: ", new_apogee_alt)
                        println("Debris seen: ", debris_removed[i,2], ", Debris removed: ", debris_removed[i,1])

                        t_last_pulse = t
                        @inbounds debris_vis_prev[i] = true
                        break # After laser was used, skip processing the other objects in this time step
                    else
                        @inbounds debris_vis_prev[i] = true
                    end
                else
                    @inbounds debris_vis_prev[i] = false
                end
            else
                @inbounds debris_vis_prev[i] = false
            end
        end

        t += dt
        push!(percentages, debris_counter / tot_debris_n)

        if mod(round(t), 50) == 0
            println("--------------------------------------------")
            println("t = ", round((t - t0) / (24 * 3600), digits=2), " days")
            println("Hit: ", count(debris_removed[:,2]))
            println("Removed: ", debris_counter)
            println(round(debris_counter / tot_debris_n * 100, digits=2), '%')

            # println(count(debris_removed[:,1] .* debris_removed[:,2]))
            if plotResults
                # Determine which debris objects are occluded
                @inbounds camera_axis = normalize([cos(view_angles[1] * pi / 180), sin(view_angles[1] * pi / 180), sin(view_angles[2] * pi / 180)])
                for i in 1:tot_debris_n
                    # Compute distance of point from camera axis
                    # Resulting distance is negative if point is on the side of Earth faced away from the camera
                    @inbounds camera_axis_dot[i] = dot(camera_axis, debris_cartesian[i,:])
                end

                # Debris that is occluded by Earth, drawn to make transition to behind Earth better
                @inbounds occluded = (camera_axis_dot .< 0) .&& .!debris_removed[:,1]
                @inbounds occluded_hit = occluded .&& debris_removed[:,2]
                @inbounds non_occluded = (camera_axis_dot .> 0) .&& .!debris_removed[:,1]
                @inbounds non_occluded_hit = non_occluded .&& debris_removed[:,2]

                # Occluded debris, not hit by laser
                @inbounds plt3d = plot(debris_cartesian[occluded, 1], debris_cartesian[occluded, 2], debris_cartesian[occluded, 3],
                    seriestype=:scatter,
                    markersize=4,
                    xlim=(-8000e3, 8000e3), ylim=(-8000e3, 8000e3), zlim=(-8000e3, 8000e3),
                    title="Space Debris Detection",
                    label="Debris fragment",
                    color=:black,
                    size=(1100, 1000),
                    camera=view_angles
                )

                # Occluded debris, hit by laser
                @inbounds scatter!(debris_cartesian[occluded_hit, 1], debris_cartesian[occluded_hit, 2], debris_cartesian[occluded_hit, 3], markersize=5, color=:red, label=false)
                
                # Spacecraft
                @inbounds scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize=10, color="green", label="Spacecraft")

                # Earth
                phi = 0:pi / 50:2 * pi
                theta = 0:pi / 100:pi
                x = [R_e * cos(t) * sin(p) for t in theta, p in phi]
                y = [R_e * sin(t) * sin(p) for t in theta, p in phi]
                z = [R_e * cos(p) for t in theta, p in phi]
                plot!(x, y, z, linetype=:surface, color=:lightblue, colorbar=false, shade=true)

                # Non-occluded debris, not hit by laser
                scatter!(debris_cartesian[non_occluded, 1], debris_cartesian[non_occluded, 2], debris_cartesian[non_occluded, 3], markersize=4, color=:black, label=false)

                # Non-occluded debris, hit by laser
                scatter!(debris_cartesian[non_occluded_hit, 1], debris_cartesian[non_occluded_hit, 2], debris_cartesian[non_occluded_hit, 3], markersize=5, color=:red, label=false)

                # Spacecraft in front of Earth
                if dot(camera_axis, position_sc) > 0
                    scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize=10, color="green", label=false)
                end

                display(plt3d)
            end
        end
    end
    increased_a_percentage = increased_a_counter / debris_counter * 100
    return (ts, percentages, increased_a_percentage)
end

@time (times, perc, perc_increased_a) = run_sim(plotResults=false)

time_required = last(times)

println("For scan time equal to ", scan_time, " s and FoV of ", FoV * 180 / pi, " deg:")
println("The time required for 50% is equal to ", round(time_required / (24 * 3600), digits=3), "days.")
println("Of which ", round(perc_increased_a, digits=3), "% have an increased semi-major axis.")
p = plot(times ./ (3600 * 24), perc .* (100 * 0.61), xlabel="Time [days]", ylabel="Removal fraction [%]", label=false)
savefig(p, string(tot_debris_n) * "-DebrisRemovalTime" * "-Cd" * string(cooldown_time) * "-fov" * string(round(FoV * 180 / pi)) * "-i" * string(round(incidence_angle * 180 / pi)) * "-r" * string(range) * "-mint" * string(min_vis_time) * ".pdf")