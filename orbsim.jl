using Base:print_array
using Pkg
# Pkg.add("LoopVectorization")
# Pkg.add("CSV")
# Pkg.add("DataFrames")

using CSV
using DataFrames
using Statistics
using LinearAlgebra
using LoopVectorization

# Constants
const R_e = 6378.137e3  # [m]
const g_0 = 9.80665  # [m/s2]
const J_2 = 0.00108263  # [-]
const mu = 3.986004418e14  # [m3/s2]
const h_collision = 789e3  # [m]
const debris_n = 10000  # total n is 21627 fragments, change this number for simulation speed

const a_collision = R_e + h_collision
const t0 = 72 * 100 * 60
const dt = 50
const distance_sc = 40e3

# Spacecraft variables
const a_sc = R_e + h_collision + distance_sc
const e_sc = 0
const M_0_sc = 0


df = CSV.read("./iridium_cosmos_result.csv", DataFrame; header=1)

# Select data that is necessary and convert to matrix
df2 = filter(row -> row.Name .== "Kosmos 2251-Collision-Fragment", df)
df3 = filter(row -> row.d_eq .< 0.1, df2)
df4 = filter(row -> row.e .< 1.0, df3)
df5 = select(df4, ["a", "e", "i", "long_asc", "arg_peri", "mean_anom"])
debris_kepler = Matrix(df5)

# Cut data set down to set number of fragments
tot_debris_n = min(debris_n, length(debris_kepler[:,1]))
println(tot_debris_n)
debris_kepler = debris_kepler[1:tot_debris_n,:]
debris_carthesian = Matrix{Float64}(undef, tot_debris_n, 3)
debris_removed = zeros(Bool, tot_debris_n)

# display(debris_kepler)

@inline function true_anom(a, e, t, M_0)
    n = sqrt(mu / a^3)
    M = n * t - M_0

    # Initial guess
    E = 0 

    # Apply newton method 3x
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

    i_sc = mean(debris_kepler[:, 3])
    RAAN_sc = mean(debris_kepler[:, 4])
    ts = Vector{Float64}(undef, 0)
    percentages = Vector{Float64}(undef, 0)
    position_sc = zeros(3)

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

    while (debris_counter / debris_n < 0.5) && (t - t0 < 7 * 24 * 3600) # Limited to 7 days for testing
        push!(ts, t)

        # Update RAAN and w due to J_2 (sc)
        RAAN_sc += RAAN_drift_sc
        w_sc += w_drift_sc

        # Update RAAN and w due to J_2 (debris)
        @inbounds debris_kepler[:, 4] += RAAN_drift
        @inbounds debris_kepler[:, 5] += w_drift

        # Compute spacecraft position
        true_anomaly_sc = true_anom(a_sc, e_sc, t, M_0_sc)
        kepler_to_cartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)
        
        # Update space debris position
        @tturbo for i = 1:tot_debris_n
            @inbounds true_anomaly_debris = true_anom(debris_kepler[i, 1], debris_kepler[i, 2], t, debris_kepler[i, 6])
            # @inbounds kepler_to_cartesian(debris_kepler[i, 1], debris_kepler[i, 2], debris_kepler[i, 5], true_anomaly_debris, debris_kepler[i, 3], debris_kepler[i, 4], view(debris_carthesian, i, :))
            

            # Convert a position in the Keplerian system to a cartesian system
            @inbounds a = debris_kepler[i, 1]
            @inbounds e = debris_kepler[i, 2]
            @inbounds w = debris_kepler[i, 5]
            @inbounds inc = debris_kepler[i, 3]
            @inbounds RAAN = debris_kepler[i, 4]

            p = a * (1 - e * e)
            r = p / (1 + e * cos(true_anomaly_debris)) # radius

            # Compute the Cartesian position vector
            @inbounds debris_carthesian[i,1] = r * (cos(RAAN) * cos(w + true_anomaly_debris) - sin(RAAN) * sin(w + true_anomaly_debris) * cos(inc))
            @inbounds debris_carthesian[i,2] = r * (sin(RAAN) * cos(w + true_anomaly_debris) + cos(RAAN) * sin(w + true_anomaly_debris) * cos(inc))
            @inbounds debris_carthesian[i,3] = r * (sin(inc) * sin(w + true_anomaly_debris))
        end

        # This is separate from the above loop because @tturbo uses vector intrinsics, which are not available for more complex functions like sqrt()
        Threads.@threads for i = 1:tot_debris_n
            @inbounds abs_distance = norm(debris_carthesian[i,:] - position_sc)
            # println(abs_distance)
            if abs_distance < 100e3
                @inbounds debris_removed[i] = true
            end
        end

        t += dt
        debris_counter = count(p -> (p .== true), debris_removed)
        push!(percentages, debris_counter / debris_n)

        println("--------------------------------------------")
        println(round((t - t0) / 3600, digits=2))
        if mod(round(t), 2) == 0
            println(debris_counter)
            println(round(debris_counter / debris_n * 100, digits=2), '%')
        end
    end
end

@time run_sim()