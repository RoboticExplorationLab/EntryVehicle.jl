#test COM on axis

using LinearAlgebra
using Plots
using ApproxFun
using Libdl

include("aero.jl")
include("dyna.jl")
include("quaternions.jl")
include("traj_plots.jl")
include("interpolation.jl")

a=1

function rk4(f, y_0, p, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(t, y_star, p)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, p)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, p)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, p)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y'
end

δ =70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.015; 0.0; 0.0]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

α = 0.0:1.0:181.0
tableFX = table_CF[:,1]
tableFY = table_CF[:,2]
tableFZ = table_CF[:,3]
tableτX = table_Cτ[:,1]
tableτY = table_Cτ[:,2]
tableτZ = table_Cτ[:,3]
order = 14
C_FX = compute_chebyshev_coefficients_aerodynamics(α, tableFX[1:182], order)
C_FY = compute_chebyshev_coefficients_aerodynamics(α, tableFY[1:182], order)
C_FZ = compute_chebyshev_coefficients_aerodynamics(α, tableFZ[1:182], order)
C_τX = compute_chebyshev_coefficients_aerodynamics(α, tableτX[1:182], order)
C_τY = compute_chebyshev_coefficients_aerodynamics(α, tableτY[1:182], order)
C_τZ = compute_chebyshev_coefficients_aerodynamics(α, tableτZ[1:182], order)



θ = 91*pi/180 #rotation angle about z body axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)

#Case when the vehicle starts pointing the surface of the planet (no initial velocity)
M = [0.0 0.0 -1.0;
    -1.0 0.0 0.0;
    0.0 1.0 0.0]
Q = mat2quat(M)
Q = qconj(Q)
v_eci = [-1.0; 5.0; 0.0]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β); -sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q)

#state ini
v_eci = [0.001; 1.0; 0.0001]*1e3
x0 = [(3389.5+125)*1e3; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]
t_sim4, Z4 = rk4(dyna_coeffoff_COM_on_axis, x0, [0.0], 0.01, [0.0, 265.0])

plot_traj(Z4)
plot_altitude(Z4, t_sim4)
plot_quaternions(Z4)
plot_attack_angle(Z4, t_sim4)

plot_ang_vel(Z4, t_sim4)

plot_vel(Z4, t_sim4)
