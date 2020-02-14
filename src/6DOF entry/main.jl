# INFORMATION
# 6DOF Entry Vehicle Model simulation
# see "aero.jl" for computation of aerodynamics coefficients
# see "dyna.jl" for the actual dynamical function 1st function called
#dyna_coeffoff_COM_on_axis(t, x, u)
# No Control in the loop here
# Sphere cone geometry here

using LinearAlgebra
using Plots
using ApproxFun
using Libdl

include("aero.jl")
include("dyna.jl")
include("quaternions.jl")
include("traj_plots.jl")
include("interpolation.jl")
include("space_mechanics.jl")


#Regular RK4 integrator ########################################################

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

################################################################################

# Vehicle geometry, aero coefficients computation and interpolation#############

# Phoenix type geometry
δ =70*pi/180 #apex angle sphere cone
r_min = 0.2 #geometry
r_cone = 1.325 #geometry
r_G = [0.001; 0.0; -0.189] #center of mass location in body frame. origin at conic base
table_CF, table_Cτ, table_damping = table_aero_spherecone(δ, r_min, r_cone, r_G)

#Sequence for interpolation of aerodynamics coefficients
α = 0.0:1.0:181.0
tableFX = table_CF[:,1]
tableFY = table_CF[:,2]
tableFZ = table_CF[:,3]
tableτX = table_Cτ[:,1]
tableτY = table_Cτ[:,2]
tableτZ = table_Cτ[:,3]
tabledx = table_damping[:, 1]
tabledy = table_damping[:, 2]
tabledz = table_damping[:, 3]
order = 14
C_FX = compute_chebyshev_coefficients_aerodynamics(α, tableFX[1:182], order)
C_FY = compute_chebyshev_coefficients_aerodynamics(α, tableFY[1:182], order)
C_FZ = compute_chebyshev_coefficients_aerodynamics(α, tableFZ[1:182], order)
C_τX = compute_chebyshev_coefficients_aerodynamics(α, tableτX[1:182], order)
C_τY = compute_chebyshev_coefficients_aerodynamics(α, tableτY[1:182], order)
C_τZ = compute_chebyshev_coefficients_aerodynamics(α, tableτZ[1:182], order)
DX = compute_chebyshev_coefficients_aerodynamics(α, tabledx[1:182], order)
DY = compute_chebyshev_coefficients_aerodynamics(α, tabledy[1:182], order)
DZ = compute_chebyshev_coefficients_aerodynamics(α, tabledz[1:182], order)

################################################################################

# Initial conditions ###########################################################

#case with non zero flight path angle
v_eci = [-1.6; 6.8; 0.0001]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β);-sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q) #quaternions initial

# Assume quasi-equatorial trajectory
# starting at 125 km
x0 = [(3389.5+125)*1e3; 0.0; 50.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]

# Integrating using RK4
t_sim4, Z4 = rk4(dyna_coeffoff_COM_on_axis, x0, [0.0], 0.01, [0.0, 180.0])

# Plotting Results #############################################################

plot_traj(Z4) #traj in X, Y plane

plot_altitude(Z4, t_sim4)

plot_quaternions(Z4, t_sim4)

#attack angle seem to converge when reaching lower level of atmosphere
plot_total_attack_angle(Z4, t_sim4)

plot_entry_profile(Z4, t_sim4)

#oscillation and stabilization?
plot_ang_vel(Z4, t_sim4)

plot_vel(Z4, t_sim4)

plot_mach_number(Z4, t_sim4)

plot_mach_number_altitude(Z4, t_sim4)

plot_specific_energy(Z4, t_sim4)
