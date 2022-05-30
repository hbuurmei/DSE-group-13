using Base:print_array
using Pkg
# Pkg.add("LoopVectorization")
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Plots")

using CSV
using DataFrames
using Statistics
using LinearAlgebra
using LoopVectorization
using Plots

# Constants
const R_e = 6378.137e3  # [m]
const g_0 = 9.80665  # [m/s2]
const J_2 = 0.00108263  # [-]
const mu = 3.986004418e14  # [m3/s2]
const h_collision = 789e3  # [m]
const debris_n = 100  # number of fragments, change this number for simulation speed

const a_collision = R_e + h_collision
const t0 = 72 * 100 * 60
const dt = 6
const distance_sc = 40e3
const target_fraction = 0.5

# Spacecraft variables
const a_sc = R_e + h_collision + distance_sc
const e_sc = 0
const M_0_sc = -0.045 * pi


df = CSV.read("./iridium_cosmos_result.csv", DataFrame; header=1)

# Select data that is necessary and convert to matrix
df = filter(row -> row.Name .== "Kosmos 2251-Collision-Fragment", df)
df = filter(row -> row.d_eq .< 0.1, df)
df = filter(row -> row.e .< 1.0, df)
debris_kepler = Matrix(select(df, ["a", "e", "i", "long_asc", "arg_peri", "mean_anom", "ID"])) # ID is used as an additional column to store true anomaly
debris_dims = Matrix(select(df, ["M"]))
debris_kepler = debris_kepler[32:33,:]
debris_dims = debris_dims[32:33,:]

# Cut data set down to set number of fragments
tot_debris_n = min(debris_n, length(debris_kepler[:,1]))
println("Number of debris objects: ", tot_debris_n)
debris_kepler = debris_kepler[1:tot_debris_n,:]
debris_true_anoms = zeros(tot_debris_n) # Include this in debris_kepler in the future
debris_cartesian = Matrix{Float64}(undef, tot_debris_n, 3)
debris_cartesian_vel = Matrix{Float64}(undef, tot_debris_n, 3)
debris_removed = zeros(Bool, tot_debris_n, 2)


@inline function calc_true_anomaly(a, e, t, M_0)
    n = sqrt(mu / a^3)
    M = n * t - M_0

    # Initial guess
    E = 0 

    # Apply newton method 10x
    for i = 1:10
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E))
    end

    # Final equation
    true_anomaly = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    return true_anomaly
end

@inline function calc_true_anomaly2(a, e, M)
    # Initial guess
    E = 0

    # Apply newton method 10x
    for i = 1:10
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E))
    end

    # Final equation
    true_anomaly = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    return true_anomaly
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

@inline function calc_vel(a, e, w, true_anomaly, i, RAAN, position)
    # Get the velocity in cartesian coordinates
    p = a * (1 - e * e)
    r = p / (1 + e * cos(true_anomaly)) # radius
    h = sqrt(mu * p)

    V_X = (position[1] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (cos(RAAN) * sin(w + true_anomaly) + sin(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Y = (position[2] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (sin(RAAN) * sin(w + true_anomaly) - cos(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Z = (position[3] * h * e / (r * p)) * sin(true_anomaly) + (h / r) * (cos(w + true_anomaly) * sin(i))

    return [V_X, V_Y, V_Z]
end

@inline function thrust_alter_orbit(debris_kepler, debris_cartesian, debris_cartesian_vel, debris_dims, thrust_dir, thrust_energy, i)
    # Compute product of a and thrust_dt
    # Based on kinetic energy and v2 = a * dt + v1
    @inbounds v1 = norm(debris_cartesian_vel[i,:])
    @inbounds dir_dv = normalize(thrust_dir) .* (sqrt(v1 * v1 + 2 * thrust_energy / debris_dims[i,1]) - v1)

    # Establish RTO (Radial, Transverse, Out-of-plane) axes (unit vectors)
    @inbounds R = normalize(debris_cartesian[i,:])
    @inbounds O = normalize(cross(R, debris_cartesian_vel[i,:]))
    T = cross(O, R)

    # Projection into RTO (Radial, Transverse, Out-of-plane)
    dir_dv_rto = zeros(3)
    dir_dv_rto[1] = dot(dir_dv, R)
    dir_dv_rto[2] = dot(dir_dv, T)
    dir_dv_rto[3] = dot(dir_dv, O)

    @inbounds sqramu = sqrt(debris_kepler[i, 1] / mu)
    @inbounds sub1e2 = 1 - debris_kepler[i, 2] * debris_kepler[i, 2]
    sqr1e2 = sqrt(sub1e2)
    @inbounds ecosf1 = debris_kepler[i, 2] * cos(debris_kepler[i, 7]) + 1
    @inbounds n = sqrt(mu / debris_kepler[i, 1]^3)

    # Gaussian perturbation formulae
    @inbounds debris_kepler[i, 1] += sqramu * 2 * debris_kepler[i, 1] / sqr1e2 * (debris_kepler[i, 2] * sin(debris_kepler[i, 7]) * dir_dv[1] + ecosf1 * dir_dv[2])
    @inbounds debris_kepler[i, 2] += sqramu * sqr1e2 * (sin(debris_kepler[i, 7]) * dir_dv[1] + (debris_kepler[i, 2] + 2 * cos(debris_kepler[i, 7]) + debris_kepler[i, 2] * cos(debris_kepler[i, 7])^2) / ecosf1 * dir_dv[2])
    @inbounds debris_kepler[i, 3] += sqramu * sqr1e2 / ecosf1 * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) * dir_dv[3]
    dRAAN = sqramu * sqr1e2 / ecosf1 * sin(debris_kepler[i, 5] + debris_kepler[i, 7]) / sin(i) * dir_dv[3]
    @inbounds debris_kepler[i, 4] += dRAAN
    @inbounds debris_kepler[i, 5] += sqramu * sqr1e2 / debris_kepler[i, 2] * (- cos(debris_kepler[i, 7]) * dir_dv[1] + (ecosf1 + 1) / ecosf1 * sin(debris_kepler[i, 7]) * dir_dv[2]) - cos(i) * dRAAN
    # @inbounds debris_kepler[i, 6] += n + sub1e2 / (n * debris_kepler[i, 1] * debris_kepler[i, 2]) * ((cos(debris_kepler[i, 7]) - 2 * debris_kepler[i, 2] / ecosf1) * dir_dv[1] - (ecosf1 + 1) / ecosf1 * sin(debris_kepler[i, 7]) * dir_dv[2])
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

function run_sim()
    debris_counter = 0
    t = t0
    w_sc = 0

    i_sc = debris_kepler[1, 3] # mean(debris_kepler[:, 3])
    RAAN_sc = debris_kepler[1, 4] # mean(debris_kepler[:, 4])
    ts = Vector{Float64}(undef, 0)
    percentages = Vector{Float64}(undef, 0)
    position_sc = zeros(3)
    debris_vis = zeros(tot_debris_n, 2) # Col1: Tot iterations visible, Col2: Number of total passes
    debris_vis_prev = zeros(Bool, tot_debris_n) # Col1: Visible prev it
    vel_sc = zeros(3)

    sizehint!(ts, 100000);
    sizehint!(percentages, 100000);

    # J_2 effect sc
    RAAN_drift_sc = J_2_RAAN(a_sc, e_sc, i_sc) * dt
    w_drift_sc = J_2_w(a_sc, e_sc, i_sc) * dt

    # J_2 effect debris
    RAAN_drift = Vector{Float64}(undef, tot_debris_n)
    w_drift = Vector{Float64}(undef, tot_debris_n)

    for i in eachindex(RAAN_drift)
        @inbounds RAAN_drift[i] = J_2_RAAN(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * dt
        @inbounds w_drift[i] = J_2_w(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * dt
    end

    # Initial mean anomalies, propate orbit until t0
    @inbounds n = sqrt.(mu / debris_kepler[:, 1].^3)
    @inbounds debris_kepler[:, 6] = transpose(t .* n) .+ debris_kepler[:, 6]

    display(debris_kepler)

    while (debris_counter / tot_debris_n < target_fraction)
        push!(ts, t - t0)

        # Update RAAN and w due to J_2 (sc)
        RAAN_sc += RAAN_drift_sc
        w_sc += w_drift_sc

        # Compute spacecraft position
        true_anomaly_sc = calc_true_anomaly(a_sc, e_sc, t, M_0_sc)
        kepler_to_cartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)

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

            # Update true anomalies
            @inbounds debris_kepler[i, 7] = calc_true_anomaly2(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 6])

            @inbounds p = debris_kepler[i, 1] * (1 - debris_kepler[i, 2] * debris_kepler[i, 2])
            @inbounds r = p / (1 + debris_kepler[i, 2] * cos(debris_kepler[i, 7])) # radius

            # Compute the Cartesian position vector
            @inbounds debris_cartesian[i, 1] = r * (cos(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) - sin(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_true_anoms[i]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 2] = r * (sin(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + debris_kepler[i, 7]) + cos(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + debris_true_anoms[i]) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 3] = r * (sin(debris_kepler[i, 3]) * sin(debris_kepler[i, 5] + debris_kepler[i, 7]))
        end

        # This is separate from the above loop because @tturbo uses vector intrinsics, which are not available for more complex functions like sqrt()
        Threads.@threads for i = 1:tot_debris_n
            @inbounds rel_pos = position_sc - debris_cartesian[i,:] # Vector from debris to spacecraft
            @inbounds abs_distance = norm(debris_cartesian[i,:] - position_sc)
            # println(abs_distance)
            if abs_distance < 100e3 # (abs_distance > 100e3) && (abs_distance < 500e3) # abs_distance < 100e3
                # Update spacecraft velocity
                @inbounds debris_velocity = calc_vel(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 5], debris_kepler[i, 7], debris_kepler[i, 3], debris_kepler[i, 4], debris_cartesian[i,:])

                # Check angle between debris tranjectory and spacecraft relative to debris
                println(sum(debris_velocity .* rel_pos) / (norm(debris_velocity) * norm(rel_pos)))
                if sum(debris_velocity .* rel_pos) / (norm(debris_velocity) * norm(rel_pos)) > 0.5 # sqrt(3) / 2
                    # Inside sphere and cone
                    @inbounds debris_vis[i,1] += 1
                    @inbounds debris_vis[i,2] += !debris_vis_prev[i] # Add only a pass if object was not detected in the previous timestep
                    @inbounds debris_vis_prev[i] = true

                    @inbounds debris_removed[i,2] = true

                    thrust_dir = - normalize(rel_pos) # Thrust in laser directions
                    thrust_energy = 100 # J
                    thrust_alter_orbit(debris_kepler, debris_cartesian, debris_cartesian_vel, debris_dims, thrust_dir, thrust_energy, i)
                else
                    # Inside sphere, but not cone
                    # println("Ineffective geometry")
                    @inbounds debris_vis_prev[i] = false
                end
            else
                # Not within range
                @inbounds debris_vis_prev[i] = false
            end
            if norm(debris_cartesian[i,:]) < (R_e + 200e3)
                @inbounds debris_removed[i,1] = true
            end
            # Update mean anomaly
            @inbounds n = sqrt.(mu / debris_kepler[i, 1].^3)
            @inbounds debris_kepler[i, 6] += n * dt
        end

        t += dt
        debris_counter = count(p -> (p .== true), debris_removed[:,2])
        push!(percentages, debris_counter / debris_n)

        if mod(round(t), 20) == 0
            println("--------------------------------------------")
            println(round((t - t0) / 3600, digits=2))
            println(debris_counter)
            println(round(debris_counter / debris_n * 100, digits=2), '%')
            display(debris_removed)

            # Visualise the orbits
            colors = [RGB(removed, 0, 0) for removed in debris_removed[:,2]]
            plt3d = plot(debris_cartesian[:, 1], debris_cartesian[:, 2], debris_cartesian[:, 3],
                seriestype=:scatter,
                markersize=4,
                xlim=(-8000e3, 8000e3), ylim=(-8000e3, 8000e3), zlim=(-8000e3, 8000e3),
                title="Space Debris Detection",
                label="Debris fragment",
                color=colors, # :black,
                size=(1100, 1000)
            )
            # Spacecraft
            scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize=5, color="green", label="Spacecraft")
            # Earth
            phi = 0:pi / 25:2 * pi
            theta = 0:pi / 50:pi
            x = [R_e * cos(t) * sin(p) for t in theta, p in phi]
            y = [R_e * sin(t) * sin(p) for t in theta, p in phi]
            z = [R_e * cos(p) for t in theta, p in phi]
            plot!(x, y, z, linetype=:surface, color=:blue, colorbar=false)

            display(plt3d)
        end
    end
    return (ts, percentages, debris_vis)
end

@time (times, perc, debris_vis_stats) = run_sim()

p = plot(times ./ (3600 * 24), perc .* (100 * 0.61), xlabel="Time [days]", ylabel="Removal fraction [%]")
savefig(p, "DebrisRemovalTime.pdf")

avg_vis_times = debris_vis_stats[:,1] .* dt ./ debris_vis_stats[:,2]
println("Average time visible: ", mean(filter(!isnan, avg_vis_times)), "s")
println("Number below 29 s: ", count(p -> (p .< 29), avg_vis_times))
h1 = histogram(avg_vis_times, xlabel="Average visibility time per pass", ylabel="Amount of debris objects", bins=40, legend=false)
savefig(h1, "DebrisVisibilityTime.pdf")
avg_vis_times[ isnan.(avg_vis_times)] .= -Inf
println(findmax(avg_vis_times))