#Functions for interpolation

using ApproxFun
using Libdl

struct Coeff_interp
    #coefficients for interpolation
    C_FX::Vector{Float64}
    C_FY::Vector{Float64}
    C_FZ::Vector{Float64}
    C_τX::Vector{Float64}
    C_τY::Vector{Float64}
    C_τZ::Vector{Float64}
    DX::Vector{Float64}
    DY::Vector{Float64}
    DZ::Vector{Float64}
end


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

function tables(geometry::SphereconeGeometry,α_min,α_max,α_step,order)
    #give interpolated coeff tables for hypersonic regime given a geometry
    #and an order in the Chebyshev polynomial
    δ = geometry.δ
    r_min = geometry.r_min
    r_cone = geometry.r_cone
    r_G = geometry.r_G
    table_CF, table_Cτ, table_damping = table_aero_spherecone(δ, r_min, r_cone, r_G)
    @show(table_CF)
    α = α_min:α_step:α_max
    tableFX = table_CF[:,1]
    tableFY = table_CF[:,2]
    tableFZ = table_CF[:,3]
    tableτX = table_Cτ[:,1]
    tableτY = table_Cτ[:,2]
    tableτZ = table_Cτ[:,3]
    tabledx = table_damping[:, 1]
    tabledy = table_damping[:, 2]
    tabledz = table_damping[:, 3]
    C_FX = compute_chebyshev_coefficients_aerodynamics(α, tableFX[1:182], order)
    C_FY = compute_chebyshev_coefficients_aerodynamics(α, tableFY[1:182], order)
    C_FZ = compute_chebyshev_coefficients_aerodynamics(α, tableFZ[1:182], order)
    C_τX = compute_chebyshev_coefficients_aerodynamics(α, tableτX[1:182], order)
    C_τY = compute_chebyshev_coefficients_aerodynamics(α, tableτY[1:182], order)
    C_τZ = compute_chebyshev_coefficients_aerodynamics(α, tableτZ[1:182], order)
    DX = compute_chebyshev_coefficients_aerodynamics(α, tabledx[1:182], order)
    DY = compute_chebyshev_coefficients_aerodynamics(α, tabledy[1:182], order)
    DZ = compute_chebyshev_coefficients_aerodynamics(α, tabledz[1:182], order)
    return C_FX, C_FY, C_FZ, C_τX, C_τY, C_τZ, DX, DY, DZ
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
