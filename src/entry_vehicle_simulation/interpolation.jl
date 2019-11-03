#Functions for interpolation

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
