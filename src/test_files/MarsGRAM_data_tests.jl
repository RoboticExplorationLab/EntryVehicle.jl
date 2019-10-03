using DataFrames
using CSV
using Plots
using ApproxFun
using LinearAlgebra
#using JLD #loading and saving julia variables, dont need that
#using ChebyshevApprox #notworking
pyplot()

#Extract MarsGRAM data after execution
#=
X = CSV.read(joinpath(pwd(), "OUTPUT.txt"))
show(X)
describe(X)
dropmissing(X)

XX = CSV.File(joinpath(pwd(), "OUTPUT.txt"))
=#
XXX = CSV.File(joinpath(pwd(),"OUTPUT.txt"); delim=' ', ignorerepeated=true)
XXX = DataFrame(XXX)

#Density with respect to Altitude for the specified conditions
plot(XXX.Denkgm3, XXX.HgtMOLA)
#Profile Temperature
plot(XXX.Temp, XXX.HgtMOLA)
plot!(XXX.Temp, XXX.HgtMOLA)

plot(XXX.Denkgm3, XXX.HgtMOLA)

#Chebyshev Polynomial Fitting using MarsGRAM data

#Simon's previous work on that (Chebyshev fitting)
function compute_chebyshev_coefficients_atmosphere_log(height, density, order)
    # We approximate this data using a Chebyshev polynomial.
    # Note that an evenly spaced grid suffers from instability for large n.
    # The easiest way around this is to use least squares with more points than coefficients,
    # instead of interpolation:
    altitude_min = height[1]
    altitude_max = height[end]
    num_nodes = size(density)[1]
    space = Chebyshev(Interval(altitude_min, altitude_max))
    # A non-default grid
    points = range(altitude_min, stop=altitude_max, length=num_nodes)
    values = log.(density)
    # Create a Vandermonde matrix by evaluating the basis at the grid
    V = zeros(Float64, num_nodes, order)
    for k = 1:order
        V[:, k] = Fun(space,[zeros(k-1);1]).(points)
    end
    atmosphere_polynomial_log = Fun(space, V\values)
    atmosphere_coefficients = atmosphere_polynomial_log.coefficients
    return atmosphere_coefficients
end

interpol_coeff = compute_chebyshev_coefficients_atmosphere_log(XXX.HgtMOLA, XXX.Denkgm3, 10)
alt_min = minimum(XXX.HgtMOLA);
alt_max = maximum(XXX.HgtMOLA);

function atmosphere_density_chebyshev(altitude)
    atmosphere_coefficients = interpol_coeff
    altitude_min = alt_min
    altitude_max = alt_max
    #atmosphere_coefficients, altitude_min, altitude_max =
    #	load_chebyshev_approximation_atmosphere_log("chebyshev_coefficients_atmosphere.jld")
    space = Chebyshev(Interval(altitude_min, altitude_max))
    density_polynomial_log = Fun(space, atmosphere_coefficients)
    density = exp(density_polynomial_log(altitude))
    return density
end

#Verification Plots

altitude = 0.0:0.5:100
dens = zeros(length(altitude))
for i = 1:length(altitude)
    dens[i] = atmosphere_density_chebyshev(altitude[i])
end
plot(dens, altitude)
ylabel!("altitude")
xlabel!("density")
