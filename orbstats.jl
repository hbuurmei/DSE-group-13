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
const view_angles = (45, 45) # Viewing angles in azimuth and altitude


# Constants
const R_e = 6378.137e3  # [m]
const g_0 = 9.80665  # [m/s2]
const J_2 = 0.00108263  # [-]
const mu = 3.986004418e14  # [m3/s2]
const h_collision = 789e3  # [m]
const debris_n = 100000  # number of fragments, change this number for simulation speed

const a_collision = R_e + h_collision
const t0 = 5 * 24 * 60 * 60 # 5 days after collision
const t_end = t0 + 20 * 24 * 60 * 60 # Run for 100 days
const dt = 6
const distance_sc = 40e3
const target_fraction = 0.5
const max_dv = 1 # Maximum dV use in gaussian perturbation equations

# Spacecraft variables
const a_sc = R_e + h_collision + distance_sc
const e_sc = 0
const M_0_sc = -0.045 * pi



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
@inline function kepler_to_cartesian(a, e, w, true_anomaly, i, RAAN, position)
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

@inline function J_2_RAAN(a, e, i)
    n = sqrt(mu / a^3)
    RAAN_dot = -1.5 * n * R_e * R_e * J_2 * cos(i) / (a * a) / (1 - e * e)^2
    return RAAN_dot
end

@inline function J_2_w(a, e, i)
    n = sqrt(mu / a^3)
    w_dot = 0.75 * n * R_e * R_e * J_2 * (4 - 5 * (sin(i))^2) / (a * a) / (1 - e * e)^2
    return w_dot
end

function run_sim(;plotResults=true)
    t = t0
    i_sc = mean(debris_kepler[:, 3]) # mean(debris_kepler[:, 3])
    w_sc = J_2_w(a_sc, e_sc, i_sc) * t0
    RAAN_sc = mean(debris_kepler[:, 4]) + J_2_RAAN(a_sc, e_sc, i_sc) * t0
    position_sc = zeros(3)
    debris_vis = zeros(tot_debris_n, 2) # Col1: Tot iterations visible, Col2: Number of total passes
    debris_vis_prev = zeros(Bool, tot_debris_n) # Col1: Visible in previous iteration
    vel_sc = zeros(3)
    camera_axis_dot = zeros(tot_debris_n)

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

    while (t < t_end)
        # Update RAAN and w due to J_2 (sc)
        RAAN_sc += RAAN_drift_sc
        w_sc += w_drift_sc

        n_sc = sqrt(mu / a_sc^3)
        true_anomaly_sc = calc_true_anomaly(a_sc, e_sc, n_sc * t + M_0_sc)
        kepler_to_cartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)

        # Update space debris position
        @tturbo for i = 1:tot_debris_n
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
        
            @inbounds rel_pos_x = position_sc[1] - debris_cartesian[i,1]
            @inbounds rel_pos_y = position_sc[2] - debris_cartesian[i,2]
            @inbounds rel_pos_z = position_sc[3] - debris_cartesian[i,3]
            @inbounds abs_distance = sqrt(rel_pos_x^2 + rel_pos_y^2 + rel_pos_z^2)

            # Get the velocity in cartesian coordinates
            @inbounds p = debris_kepler[i, 1] * (1 - debris_kepler[i, 2] * debris_kepler[i, 2])
            @inbounds r = p / (1 + debris_kepler[i, 2] * cos(debris_kepler[i, 7])) # radius
            h = sqrt(mu * p)

            @inbounds debris_cartesian_vel[i,1] = (debris_cartesian[i,1] * h * debris_kepler[i, 2] / (r * p)) * sin(debris_kepler[i, 7]) - (h / r) * (cos(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) + sin(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian_vel[i,2] = (debris_cartesian[i,2] * h * debris_kepler[i, 2] / (r * p)) * sin(debris_kepler[i, 7]) - (h / r) * (sin(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) - cos(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian_vel[i,3] = (debris_cartesian[i,3] * h * debris_kepler[i, 2] / (r * p)) * sin(debris_kepler[i, 7]) + (h / r) * (cos(debris_kepler[i, 5] + debris_kepler[i, 7]) * sin(debris_kepler[i, 3]))

            in_range = (100e3 < abs_distance < 500e3)
            @inbounds vel_norm = sqrt(debris_cartesian_vel[i,1]^2 + debris_cartesian_vel[i,2]^2 + debris_cartesian_vel[i,3]^2)
            @inbounds rel_pos_vel_pos_dot = debris_cartesian_vel[i,1] * rel_pos_x + debris_cartesian_vel[i,2] * rel_pos_y + debris_cartesian_vel[i,3] * rel_pos_z
            in_angle = (rel_pos_vel_pos_dot / (vel_norm * abs_distance) > (sqrt(3) / 2))

            @inbounds debris_vis[i,1] += in_range * in_angle
            @inbounds debris_vis[i,2] += (debris_vis_prev[i] ? false : true) * in_range * in_angle
            @inbounds debris_vis_prev[i] = in_range * in_angle
        end

        t += dt

        if mod(round(t), 50) == 0
            println("t = ", round((t - t0) / (24 * 3600), digits=2), " days")

            if plotResults
                # Determine which debris objects are occluded
                camera_axis = normalize([cos(view_angles[1] * pi / 180), sin(view_angles[1] * pi / 180), sin(view_angles[2] * pi / 180)])
                for i in 1:tot_debris_n
                    # Compute distance of point from camera axis
                    # Resulting distance is negative if point is on the side of Earth faced away from the camera
                    camera_axis_dot[i] = dot(camera_axis, debris_cartesian[i,:])
                end

                # Debris that is occluded by Earth, drawn to make transition to behind Earth better
                occluded = (camera_axis_dot .< 0)
                non_occluded = (camera_axis_dot .> 0)
                
                # Occluded debris
                plt3d = plot(debris_cartesian[.!occluded, 1], debris_cartesian[.!occluded, 2], debris_cartesian[.!occluded, 3],
                    seriestype=:scatter,
                    markersize=4,
                    xlim=(-8000e3, 8000e3), ylim=(-8000e3, 8000e3), zlim=(-8000e3, 8000e3),
                    title="Space Debris Detection",
                    label="Debris fragment",
                    color=:black,
                    size=(1100, 1000),
                    camera=view_angles
                )

                # Spacecraft
                scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize=10, color="green", label="Spacecraft")
                # Earth
                phi = 0:pi / 50:2 * pi
                theta = 0:pi / 100:pi
                x = [R_e * cos(t) * sin(p) for t in theta, p in phi]
                y = [R_e * sin(t) * sin(p) for t in theta, p in phi]
                z = [R_e * cos(p) for t in theta, p in phi]
                plot!(x, y, z, linetype=:surface, color=:lightblue, colorbar=false, shade=true)

                # Debris in front of Earth
                scatter!(debris_cartesian[.!non_occluded, 1], debris_cartesian[.!non_occluded, 2], debris_cartesian[.!non_occluded, 3], markersize=4, color=:black, label=false)

                # Spacecraft in front of Earth
                if dot(camera_axis, position_sc) > 0
                    scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize=10, color="green", label=false)
                end

                display(plt3d)
            end
        end
    end
    return (debris_vis)
end

@time (debris_vis_stats) = run_sim(plotResults=false)

avg_vis_times = debris_vis_stats[:,1] .* dt ./ debris_vis_stats[:,2]
println("Average time visible: ", mean(filter(!isnan, avg_vis_times)), "s")
println("% of particles with visibility time below 30 s: ", count(p -> (p .< 30), avg_vis_times) / tot_debris_n * 100)
h1 = histogram(filter(vis_time -> vis_time < 20, avg_vis_times), xlabel="Average visibility time per pass", ylabel="Amount of debris objects", bins=40, legend=false)
savefig(h1, "DebrisVisibilityTime.pdf")