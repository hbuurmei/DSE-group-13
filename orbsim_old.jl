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
const debris_n = 1000  # number of fragments, change this number for simulation speed

const a_collision = R_e + h_collision
const t0 = 72 * 100 * 60.
const dt = 6
const distance_sc = 30e3  # Orbital altitude of spacecraft above collision altitude
const target_fraction = 0.5

# Spacecraft variables
const a_sc = R_e + h_collision + distance_sc
const e_sc = 0
const M_0_sc = 0


df = CSV.read("./iridium_cosmos_result.csv", DataFrame; header=1)

# Select data that is necessary and convert to matrix
df = filter(row -> row.Name .== "Kosmos 2251-Collision-Fragment", df)
df = filter(row -> row.d_eq .< 0.1, df)  # Select debris smaller than 10 cm (1 cm already set as a minimum)
df = filter(row -> row.e .< 1.0, df)  # Select debris which does not have a hyperbolic trajectory
df = filter(row -> row.a * (1 - row.e) .> (R_e + 200e3), df) # Filter out all that already have a low enough perigee
# df = filter(row -> row.a .> (a_collision - 60e3), df)
df = select(df, ["a", "e", "i", "long_asc", "arg_peri", "mean_anom"])
debris_kepler = Matrix(df)

# Cut data set down to set number of fragments
tot_debris_n = min(debris_n, length(debris_kepler[:,1]))
println("# debris considered: ", tot_debris_n)
debris_kepler = debris_kepler[1:tot_debris_n,:]
debris_cartesian = Matrix{Float64}(undef, tot_debris_n, 3)
debris_removed = zeros(Bool, tot_debris_n)


@inline function true_anom(a, e, t, M_0)
    n = sqrt(mu / a^3)
    M = n * t + M_0  # Mean anomaly. t is the time after same epoch as M_0, so don't subtract t0!

    # Initial guess
    E = 0.

    # Apply newton method 10x to find eccentric anomaly
    for i = 1:10
        E = E - (E - e * sin(E) - M) / (1 - e * cos(E))
    end

    # Final equation for true anomaly
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

@inline function getVelocity(a, e, w, true_anomaly, i, RAAN, position, velocity)
    # Get the velocity in cartesian coordinates
    p = a * (1 - e * e)
    r = p / (1 + e * cos(true_anomaly)) # radius
    h = sqrt(mu * p)

    V_X = (position[1] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (cos(RAAN) * sin(w + true_anomaly) + sin(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Y = (position[2] * h * e / (r * p)) * sin(true_anomaly) - (h / r) * (sin(RAAN) * sin(w + true_anomaly) - cos(RAAN) * cos(w + true_anomaly) * cos(i))
    V_Z = (position[3] * h * e / (r * p)) * sin(true_anomaly) + (h / r) * (cos(w + true_anomaly) * sin(i))

    velocity[1] = V_X
    velocity[2] = V_Y
    velocity[3] = V_Z
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
    debris_counter = 0
    increased_a_counter = 0
    t = t0

    i_sc = mean(debris_kepler[:, 3])
    w_sc = J_2_w(a_sc, e_sc, i_sc) * t
    RAAN_sc = mean(debris_kepler[:, 4]) + J_2_RAAN(a_sc, e_sc, i_sc) * t
    ts = Vector{Float64}(undef, 0)
    percentages = Vector{Float64}(undef, 0)
    position_sc = zeros(3)
    vel_sc = zeros(3)
    vel_debris = zeros(3)

    sizehint!(ts, 100000);
    sizehint!(percentages, 100000);

    # Take into account the J2 effect on the initial time
    for i = 1:tot_debris_n
        debris_kepler[i, 4] += J_2_RAAN(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * t
        debris_kepler[i, 5] += J_2_w(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 3]) * t
    end

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

    while (debris_counter / tot_debris_n < target_fraction)
        push!(ts, t - t0)

        # Update RAAN and w due to J_2 (sc)
        RAAN_sc += RAAN_drift_sc
        w_sc += w_drift_sc

        # Compute spacecraft position
        true_anomaly_sc = true_anom(a_sc, e_sc, t, M_0_sc)
        kepler_to_cartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)

        # Update space debris position
        @turbo for i = 1:tot_debris_n
            # Update RAAN and w due to J_2 (debris)
            @inbounds debris_kepler[i, 4] += RAAN_drift[i]
            @inbounds debris_kepler[i, 5] += w_drift[i]

            @inbounds true_anomaly_debris = true_anom(debris_kepler[i, 1], debris_kepler[i, 2], t, debris_kepler[i, 6])

            # left here for readability
            # @inbounds a = debris_kepler[i, 1]
            # @inbounds e = debris_kepler[i, 2]
            # @inbounds w = debris_kepler[i, 5]
            # @inbounds inc = debris_kepler[i, 3]
            # @inbounds RAAN = debris_kepler[i, 4]

            p = debris_kepler[i, 1] * (1 - debris_kepler[i, 2] * debris_kepler[i, 2])
            r = p / (1 + debris_kepler[i, 2] * cos(true_anomaly_debris)) # radius

            # Compute the Cartesian position vector
            @inbounds debris_cartesian[i, 1] = r * (cos(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + true_anomaly_debris) - sin(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + true_anomaly_debris) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 2] = r * (sin(debris_kepler[i, 4]) * cos(debris_kepler[i, 5] + true_anomaly_debris) + cos(debris_kepler[i, 4]) * sin(debris_kepler[i, 5] + true_anomaly_debris) * cos(debris_kepler[i, 3]))
            @inbounds debris_cartesian[i, 3] = r * (sin(debris_kepler[i, 3]) * sin(debris_kepler[i, 5] + true_anomaly_debris))
        end

        # Update spacecraft velocity
        getVelocity(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc, vel_sc)

        # This is separate from the above loop because @tturbo uses vector intrinsics, which are not available for more complex functions like sqrt()
        Threads.@threads for i = 1:tot_debris_n
            @inbounds rel_pos = position_sc - debris_cartesian[i, :]
            @inbounds abs_distance = norm(rel_pos)
            if abs_distance < 100e3  && !debris_removed[i]
                @inbounds true_anomaly_debris = true_anom(debris_kepler[i, 1], debris_kepler[i, 2], t, debris_kepler[i, 6])
                getVelocity(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 5], true_anomaly_debris,
                                         debris_kepler[i, 3], debris_kepler[i, 4], debris_cartesian[i, :], vel_debris)
                if sum(vel_debris .* rel_pos) / (norm(vel_debris) * norm(rel_pos)) > 0.5
                    @inbounds debris_removed[i] = true
                    if debris_kepler[i, 1] > a_collision
                        increased_a_counter += 1
                    end
                else
                    # println("Ineffective geometry")
                end
            end
        end

        t += dt
        debris_counter = count(p -> (p .== true), debris_removed)
        push!(percentages, debris_counter / tot_debris_n)

        if mod(round(t), dt) == 0  # Smoothest visualisation is for dt
            println("--------------------------------------------")
            println(round((t - t0) / 3600, digits=2))
            println(debris_counter)
            println(round(debris_counter / tot_debris_n * 100, digits=2), '%')

            if plotResults
                # Visualise the orbits
                plt3d = plot(debris_cartesian[:, 1], debris_cartesian[:, 2], debris_cartesian[:, 3],
                    seriestype=:scatter,
                    markersize = 3,
                    xlim = (-9000e3, 9000e3), ylim = (-9000e3, 9000e3), zlim = (-9000e3, 9000e3),
                    title = "Space Debris Detection",
                    label = "Debris fragment",
                    color = "black",
                    size = (1200, 1000),  # width, height
                    alpha = 1.0,
                    #camera = (45, 45),  # horizontal angle, vertical angle
                )

                # Spacecraft
                scatter!([position_sc[1]], [position_sc[2]], [position_sc[3]], markersize = 15, color = "green", label = "Spacecraft", alpha = 1.0)

                # Earth
                phi = 0:pi/50:2*pi
                theta = 0:pi/100:pi
                x = [R_e*cos(t)*sin(p) for t in theta, p in phi]
                y = [R_e*sin(t)*sin(p) for t in theta, p in phi]
                z = [R_e*cos(p) for t in theta, p in phi]
                display(surface!(x, y, z, color = "blue", colorbar = false, alpha=0.1))
                # display(plt3d)
            end
        end
    end
    increased_a_percentage = increased_a_counter/debris_counter*100
    return (ts, percentages, increased_a_percentage)
end

@time (times, perc, perc_increased_a) = run_sim(plotResults=true)

time_required = times[end]/(3600*24)  # In days

println("The time required for 50% is equal to ", round(time_required, digits=3), "days.")
println("Of which ", round(perc_increased_a, digits=3), "% have an increased semi-major axis.")

# Plot result
# p = plot(times./(3600*24), perc.*100, xlabel="Time [days]", ylabel="Removal fraction [%]", legend = false)
# savefig(p, "DebrisRemovalTime.pdf")
