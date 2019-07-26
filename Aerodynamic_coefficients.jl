using LinearAlgebra
using HCubature
using StaticArrays

#Coefficients computed in the frame 2 (rotated wrt body such that β = 0)
function  cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)
    #delta is the half-angle of the cone
    #alpha is the angle of attack
    u = x[1]
    v = x[2]
    sinθ = cos(α)*sin(δ)+sin(α)*cos(v)*cos(δ) #sin or cos depends on Y and Z in fact (cos means Y)
    sinθ = max(0, sinθ)
    n̂ = [-sin(δ); -cos(v)*cos(δ); sin(v)*cos(δ)]
    dA = u/sin(δ)
    Cp = 2*(sinθ)^2
    df = -Cp*dA.*n̂ #minus because the paper considered the normal vector pointing inward
    dCx = dot(df, [1;0;0])
    dCy = dot(df, [0;1;0])
    dCz = dot(df, [0;0;1])
    h1 = r_base/tan(δ)
    #h2 = r_top/tan(δ)
    h = h1 #total length of the cone
    r = [h-u/tan(δ); u*cos(v); -u*sin(v)]
    dCl = dot(cross(r, df), [1;0;0])
    dCm = dot(cross(r, df), [0;1;0])
    dCn = dot(cross(r, df), [0;0;1])
    return dCx, dCy, dCz, dCl, dCm, dCn
end

function  sphere_integrand(x::SVector{2,Float64}, α, r_n)
    #delta is the half-angle of the cone
    #alpha is the angle of attack
    u = x[1]
    v = x[2]
    sinθ = +cos(α)*u/r_n+sin(α)*cos(v)*sqrt(1-(u/r_n)^2) #check minuses or pluses here
    sinθ = max(0, sinθ)
    n̂ = [u/r_n; sqrt(1-(u/r_n)^2)*cos(v); -sqrt(1-(u/r_n)^2)*sin(v)]
    dA = r_n
    Cp = 2*(sinθ)^2
    df = -Cp*dA.*n̂ #minus because the paper considered the normal vector pointing inward
    dCx = dot(df, [1;0;0]) #X2
    dCy = dot(df, [0;1;0]) #Y2
    dCz = dot(df, [0;0;1]) #Z2
    r = [u; sqrt(r_n^2-u^2)*cos(v); -sqrt(r_n^2-u^2)*sin(v)] #check here as well then
    dCl = dot(cross(r, df), [1;0;0])
    dCm = dot(cross(r, df), [0;1;0])
    dCn = dot(cross(r, df), [0;0;1])
    return dCx, dCy, dCz, dCl, dCm, dCn
end

function cone_aero_coefficients(δ, α, r_base, r_top)

    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000

    # Integration Bounds

    r_min = r_top
    r_max = r_base
    start_bound = [r_min, 0.0]
    end_bound = [r_max, 2 * pi]

    # Define functions for Integration

    function cone_int_1(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[1]
    end
    function cone_int_2(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[2]
    end
    function cone_int_3(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[3]
    end
    function cone_int_4(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[4]
    end
    function cone_int_5(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[5]
    end
    function cone_int_6(x::SVector{2,Float64})
        return cone_integrand(x::SVector{2,Float64}, δ, α, r_base, r_top)[6]
    end

    # Integration aerodynamic force coefficients

    integral_FX, error_FX = hcubature(cone_int_1, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FX = integral_FX
    integral_FY, error_FY = hcubature(cone_int_2, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FY = integral_FY
    integral_FZ, error_FZ = hcubature(cone_int_3, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FZ = integral_FZ

    # Integration aerodynamic torque coefficients

    integral_τl, error_τl = hcubature(cone_int_4, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τl = integral_τl
    integral_τm, error_τm = hcubature(cone_int_5, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τm = integral_τm
    integral_τn, error_τn = hcubature(cone_int_6, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τn = integral_τn

    C_F = [C_FX, C_FY, C_FZ]
    C_τ = [C_τl, C_τm, C_τn]
    return C_F, C_τ
end

function sphere_aero_coefficients(α, r_n)
    #r_n = 0.5847608800326173 in our test case
    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000

    # Integration Bounds
    r_min = 0.2
    x_min = sqrt(r_n^2-r_min^2)
    x_max = r_n
    start_bound = [x_min, 0.0]
    end_bound = [x_max, 2 * pi]

    # Define functions for Integration

    function sphere_int_1(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[1]
    end
    function sphere_int_2(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[2]
    end
    function sphere_int_3(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[3]
    end
    function sphere_int_4(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[4]
    end
    function sphere_int_5(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[5]
    end
    function sphere_int_6(x::SVector{2,Float64})
        return sphere_integrand(x::SVector{2,Float64}, α, r_n)[6]
    end

    # Integration aerodynamic force coefficients

    integral_FX, error_FX = hcubature(sphere_int_1, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FX = integral_FX
    integral_FY, error_FY = hcubature(sphere_int_2, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FY = integral_FY
    integral_FZ, error_FZ = hcubature(sphere_int_3, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FZ = integral_FZ

    # Integration aerodynamic torque coefficients

    integral_τl, error_τl = hcubature(sphere_int_4, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τl = integral_τl
    integral_τm, error_τm = hcubature(sphere_int_5, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τm = integral_τm
    integral_τn, error_τn = hcubature(sphere_int_6, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τn = integral_τn

    C_F = [C_FX, C_FY, C_FZ]
    C_τ = [C_τl, C_τm, C_τn]
    return C_F, C_τ
end

function table_aero_pos(δ, r_base, r_top)
    T_X2 = []
    T_Y2 = []
    T_n2 = []
    for i=0:1:180
        C_F, C_τ = cone_aero_coefficients(δ, i*pi/180, r_base, r_top)
        append!(T_X2, C_F[1])
        append!(T_Y2, C_F[2])
        append!(T_n2, C_τ[3])
    end
    return T_X2, T_Y2, T_n2
end

function table_aero_neg(δ, r_base, r_top)
    T_X2 = []
    T_Y2 = []
    T_n2 = []
    for i=0:1:180
        C_F, C_τ = cone_aero_coefficients(δ, -i*pi/180, r_base, r_top)
        append!(T_X2, C_F[1])
        append!(T_Y2, C_F[2])
        append!(T_n2, C_τ[3])
    end
    return T_X2, T_Y2, T_n2
end

#Trying to compute with Beta as well for comparison with the paper

function  cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)
    #delta is the half-angle of the cone
    #alpha is the angle of attack
    u = x[1]
    v = x[2]
    #first is mine
    sinθ = sin(α)*cos(β)*cos(v)*cos(δ) + sin(β)*sin(v)*cos(δ) + sin(δ)*cos(α)*cos(β)
    sinθ = max(0, sinθ)
    n̂ = [cos(v)*cos(δ); sin(v)*cos(δ); sin(δ)]
    #=sinθ = cos(α)*sin(δ)*cos(β)+sin(β)*cos(v)*cos(δ)-sin(α)*cos(β)*sin(v)*cos(δ) #sin or cos depends on Y and Z in fact (cos means Y)
    sinθ = max(0, sinθ)
    n̂ = [-sin(δ); -cos(v)*cos(δ); sin(v)*cos(δ)]=#
    dA = u/sin(δ)
    Cp = 2*(sinθ)^2
    df = -Cp*dA.*n̂
    dCx = dot(df, [1;0;0])
    dCy = dot(df, [0;1;0])
    dCz = dot(df, [0;0;1])
    h1 = r_base/tan(δ)
    #h2 = r_top/tan(δ)
    h = h1 #total length of the cone
    r = [u*cos(v); u*sin(v); h-u/tan(δ)]
    dCl = dot(cross(r, df), [1;0;0])
    dCm = dot(cross(r, df), [0;1;0])
    dCn = dot(cross(r, df), [0;0;1])
    return dCx, dCy, dCz, dCl, dCm, dCn
end

function cone_aero_coefficients2(δ, α, β, r_base, r_top)

    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000

    # Integration Bounds

    r_min = r_top
    r_max = r_base
    start_bound = [r_min, 0.0]
    end_bound = [r_max, 2 * pi]

    # Define functions for Integration

    function cone_int_1(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[1]
    end
    function cone_int_2(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[2]
    end
    function cone_int_3(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[3]
    end
    function cone_int_4(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[4]
    end
    function cone_int_5(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[5]
    end
    function cone_int_6(x::SVector{2,Float64})
        return cone_integrand2(x::SVector{2,Float64}, δ, α, β, r_base, r_top)[6]
    end

    # Integration aerodynamic force coefficients

    integral_FX, error_FX = hcubature(cone_int_1, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FX = integral_FX
    integral_FY, error_FY = hcubature(cone_int_2, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FY = integral_FY
    integral_FZ, error_FZ = hcubature(cone_int_3, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_FZ = integral_FZ

    # Integration aerodynamic torque coefficients

    integral_τl, error_τl = hcubature(cone_int_4, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τl = integral_τl
    integral_τm, error_τm = hcubature(cone_int_5, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τm = integral_τm
    integral_τn, error_τn = hcubature(cone_int_6, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    C_τn = integral_τn

    C_F = [C_FX, C_FY, C_FZ]
    C_τ = [C_τl, C_τm, C_τn]
    return C_F, C_τ
end

function table_cone(δ, r_base, r_top)
    T_X = zeros(91, 91)
    T_Y = zeros(91, 91)
    T_Z = zeros(91, 91)
    T_l = zeros(91, 91)
    T_m = zeros(91, 91)
    T_n = zeros(91, 91)
    for i=0:1:90
        for j=0:1:90
        C_F, C_τ = cone_aero_coefficients2(δ, i*pi/180, j*pi/180, r_base, r_top)
        T_X[i+1, j+1] = -C_F[1]
        T_Y[i+1, j+1] = -C_F[2]
        T_Z[i+1, j+1] = -C_F[3]
        T_l[i+1, j+1] = -C_τ[1]
        T_m[i+1, j+1] = -C_τ[2]
        T_n[i+1, j+1] = -C_τ[3]
        end
    end
    return T_X, T_Y, T_Z, T_l, T_m, T_n
end

function l_ref_integrand(x::SVector{2,Float64}, δ, r_base)
    u = x[1]
    v = x[2]
    return sqrt((r_base/tan(δ)-u/tan(δ))^2+u^2)*u
end

function l_ref_cone(δ, r_base, r_top)
    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000

    # Integration Bounds

    r_min = r_top
    r_max = r_base
    start_bound = [r_min, 0.0]
    end_bound = [r_max, 2 * pi]

    function integrand(x::SVector{2,Float64})
        return l_ref_integrand(x::SVector{2,Float64}, δ, r_base)
    end

    l_ref, error_τl = hcubature(integrand, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    return l_ref

end

function check_coeff(δ, r_base, r_top)
        T_X = zeros(61)
        T_Z = zeros(61)
        T_n = zeros(61)
        for i=0:1:60
            C_F, C_τ = cone_aero_coefficients2(δ, i*pi/180, 20*pi/180, r_base, r_top)
            T_X[i+1] = C_F[1]
            T_Z[i+1] = C_F[3]
            T_n[i+1] = C_τ[3]
        end
        return T_X/5.3, T_Z/5.3, T_n
    end
