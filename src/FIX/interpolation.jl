#Functions for interpolation

#Atmopsheric density
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

function atmosphere_density_chebyshev(altitude, alt_min, alt_max)
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

#Aerodynamics Coefficients

function compute_chebyshev_coefficients_aerodynamics(α, table, order)
    # We approximate this data using a Chebyshev polynomial.
    # Note that an evenly spaced grid suffers from instability for large n.
    # The easiest way around this is to use least squares with more points than coefficients,
    # instead of interpolation:
    α_min = α[1]
    α_max = α[end]
    num_nodes = size(table)[1]
    space = Chebyshev(Interval(α_min, α_max))
    # A non-default grid
    points = range(α_min, stop=α_max, length=num_nodes)
    values = table
    # Create a Vandermonde matrix by evaluating the basis at the grid
    V = zeros(Float64, num_nodes, order)
    for k = 1:order
        V[:, k] = Fun(space,[zeros(k-1);1]).(points)
    end
    table_polynomial_coeff = Fun(space, V\values)
    return table_polynomial_coeff.coefficients
end

function table_aero_chebyshev(α, α_min, α_max, COEFF)
    polynomial_coefficients = COEFF
    min = α_min
    max = α_max
    #atmosphere_coefficients, altitude_min, altitude_max =
    #	load_chebyshev_approximation_atmosphere_log("chebyshev_coefficients_atmosphere.jld")
    space = Chebyshev(Interval(min, max))
    table = Fun(space, polynomial_coefficients)
    table_result = table(α)
    return table_result
end


#=
using ApproxFun
using Libdl


δ =70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.0; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

α = 0.0:1.0:181.0
table = table_Cτ[:,2]
order = 14
COEFF = compute_chebyshev_coefficients_aerodynamics(α, table[1:182], order)


C_chebyshev = [table_aero_chebyshev(α_e, 0.0, 181.0) for α_e in α]

Plots.plot(α, table[1:182])
Plots.plot!(α, C_chebyshev) =#
