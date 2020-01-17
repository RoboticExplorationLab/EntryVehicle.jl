#test COM on axis

using LinearAlgebra
using Plots
using ApproxFun
using Libdl
using MosekTools
using Mosek
using MathOptInterface
using JuMP

include("aero.jl")
include("aero_full_on.jl")
include("dyna.jl")
include("quaternions.jl")
include("traj_plots.jl")
include("interpolation.jl")
include("ellipsoid.jl")
include("space_mechanics.jl")


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
r_cone = 1.325
r_G = [0.001; 0.0; -0.189]
table_CF, table_Cτ, table_damping = table_aero_spherecone(δ, r_min, r_cone, r_G)

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


θ = 0.0*pi/180 #rotation angle about z body axis
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

#case with non zero flight path angle
v_eci = [-1.6; 6.8; 0.0001]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β);-sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q)

#corresponding flight-path angle
r_eci = [(3389.5+125)*1e3; 0.0; 50.0]
r_entry3DOF = ecef2entry3DOF([r_eci; v_eci])
γ = r_entry3DOF[5]*180/pi

#state ini
v_eci = [-0.001; 1.0; 0.0001]*1e3
x0 = [(3389.5+125)*1e3; 0.0; 50.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]
t_sim4, Z4 = rk4(dyna_coeffoff_COM_on_axis, x0, [0.0], 0.01, [0.0, 100.0])
t_sim4, Z4 = rk4(dyna_coeffon_COM_on_axis, x0, [0.0], 0.05, [0.0, 230.0])


plot_traj(Z4)
plot_altitude(Z4, t_sim4)
savefig("alt")
plot_quaternions(Z4, t_sim4)
savefig("quat")
#plot_attack_angle(Z4, t_sim4)
plot_total_attack_angle(Z4, t_sim4)
savefig("attack_angle")

plot_entry_profile(Z4, t_sim4)
savefig("entry_profile")
plot_ang_vel(Z4, t_sim4)
savefig("ang_vel")
plot_vel(Z4, t_sim4)
savefig("vel")
plot_mach_number(Z4, t_sim4)
plot_mach_number_altitude(Z4, t_sim4)
savefig("mach")
plot_specific_energy(Z4, t_sim4)

Plots.plot!(t_sim4, Z4[10, :])

f(u) = dyna_coeffoff_COM_on_axis(3.0, x0, [0.0])
A = ForwardDiff.jacobian(f, x0)

a=1

#############################################################
######## Uncertainty Test ###################################
#############################################################

x0_12 = [(3389.5+125)*1e3/(1e3*3389.5); 0.0; 50.0/(1e3*3389.5); 0.0; 0.0; 0.0; v_eci/(1e3*7.00); 0.0; 0.0; 0.0]#because the center with respect to himself is 0
Q0 = Diagonal([(50.0/(3389.5*1e3))^2;(50.0/(3389.5*1e3))^2;(0.1/(3389.5*1e3))^2; 0.005^2; 0.005^2; 0.005^2; (1e-3/(1e3*7.00))^2; (1e-3/(1e3*7.00))^2; (1e-4/(1e3*7.00))^2; (1e-40)^2; (1e-40)^2;(1e-40)^2]) #here we assume the change in a(vector part of the change in quaternion)

Q0 = Diagonal(ones(12)*(1e-10)^2)
Q0 = Matrix(Q0)
A0 = inv(sqrt(Q0))
b0 = -A0*x0_12

Alist, blist, centerlist, XX, ref, TT = ellipse_propagation(A0, b0, 0.0, 100.0, 1.0)

uncertainty(Alist)

function plot_traj_center(centerlist, T)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]*3389.5*1e3
        Y[j] = centerlist[2, j]*3389.5*1e3
    end
    Plots.scatter(X, Y)
end

function X_lims(X)
    n, m = size(X)
    lim = zeros(13, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
end

function XX_lims(XX)
    n, m, time = size(XX)
    lims = zeros(13, 2, time)
    for i=1:1:time
        lims[:, :, i] = X_lims(XX[:, :, i])
    end
    return lims
end

function uncertainty(Alist)
    t = length(Alist[1, 1, :])
    U = zeros(t)
    for i=1:1:t
        U[i] = tr(inv(Alist[:, :, i]))
    end
    Plots.plot(U)
    Plots.xlabel!("time")
    Plots.ylabel!("uncertainty matrix trace")
end

function euler2quat_center(centerlist, ref)
    n, t = size(centerlist)
    center_X13 = zeros(n+1, t)
    for i=1:1:t
        center_X13[:, i] = point12_13(centerlist[:, i], ref[:, i])
    end
    return center_X13
end

center13 = euler2quat_center(centerlist, ref)

Plots.scatter(TT, center13[1, :]*1e3*3389.5)

plot_traj_center(centerlist, TT)

Plots.scatter(centerlist[1, :]*3389.5*1e3, centerlist[2, :]*3389.5*1e3)

Plots.scatter(TT[1:80], centerlist[1, 1:80]*3389.5*1e3)
Plots.scatter(TT, centerlist[6, :])
Plots.scatter!(TT, centerlist[9, :]*1e3*7.00)
Plots.scatter(TT, centerlist[10, :])

Plots.plot(TT, centerlist[1, :]*3389.5*1e3)
Plots.plot!(t_sim4, Z4[1, :])


uncertainty(Alist)


lim = XX_lims(XX)

Plots.scatter(lim[1, 1, :]*3389.5*1e3)
Plots.scatter!(lim[1, 2, :]*3389.5*1e3)

Plots.scatter(lim[1, 1, :]/(1e14))
Plots.scatter!(lim[1, 2, :]/(1e14))

QQ = uncert_mat(Alist)

F = inv(Alist[:, :, end])
W = eigen(F)


################################################################################
###########################MONTE CARLO##########################################
################################################################################

using Distributions
using Random

function generate_samples(x_0, Q, M)
    #M number of samples
    n = length(x_0)
    univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    D = Product(univariate_D_vector)
    X_samples = zeros(n, M)
    rand!(D, X_samples)
    return X_samples
end

function point12_13(X)
    X13 = zeros(length(X)+1)
    X13[1:3] = X[1:3]
    X13[4:7] = euler2quat(X[4:6])
    X13[8:13] = X[7:12]
    return X13
end


function prop_MC_entry(X_samples, t_start, t_end, dt, p, model)
    n, M = size(X_samples)
    #saveAT = 1.0
    TT = t_start:dt:t_end
    traj = zeros(n+1, length(TT), M)
    for i=1:1:M
        @show(i)
        x_ini = point12_13(X_samples[:, i])
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(model, x_ini, p, dt, [t_start, t_end])
        traj[:, :, i] = Z
    end
    return traj, TT
end


θ = 20.0*pi/180.0
M = [cos(θ) -sin(θ) 0.0;
     sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0]
Q = mat2quat(M)
Q = qconj(Q)

x0 = [(3389.5+125)*1e3; 0.0; 50.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]
e = quat2euler(Q)
norm(e)
e/(norm(e))
x0_12 = [(3389.5+125)*1e3; 0.0; 50.0; e; v_eci; 0.0; 0.0; 0.0]
Q0 = Diagonal([(50.0)^2;(50.0)^2;(0.1)^2; 0.005^2; 0.005^2; 0.005^2; (1e-3)^2; (1e-3)^2; (1e-4)^2; (1e-40)^2; (1e-40)^2;(1e-40)^2])
X_samples = generate_samples2(x0_12, Q0, 500)
p = [0.0]
traj, TT = prop_MC_entry(X_samples, 0.0, 100.0, 0.01, p, dyna_coeffoff_COM_on_axis)

a = 1.0

function mean_var_MC(traj)
    n, t, M = size(traj)
    avg = zeros(n, t)
    var = zeros(n, n, t)
    for i=1:1:t
        @show(i)
        S = zeros(n)
        V = zeros(n, n)
        for j=1:1:M
            S += traj[:, i, j]
        end
        avg[:, i] = S/M
        for k =1:1:M
            V+= (traj[:, i, k]-avg[:, i])*(traj[:, i, k]-avg[:, i])'
        end
        var[:, :, i] = V/M
    end
    return avg, var
end

avg, var = mean_var_MC(traj)

Plots.scatter(0.0:0.01:150.0, avg[3, :])
Plots.plot!(t_sim4, Z4[3, :])


S = [inv(Alist[:, :, i])[1, 1] for i=1:1:length(0:1:100)]
Plots.plot(T, (centerlist[1, :]*3389.5).-3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(T, (centerlist[1, :]-S.-1.0)*3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(T, (centerlist[1, :]+S.-1.0)*3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(TT, ((traj[1, :, 76])-avg[1, :]), linewidth=0.1, legend = false)
V = [sqrt(var[1, 1, i]) for i=1:1:length(TT)]
Plots.plot!(TT, (avg[1, :]-Z4[1, :]-3*V))
Plots.plot!(TT, (avg[1, :]-Z4[1, :]+3*V))

Plots.plot(TT[1:5700], (-3*V[1:5700]))
Plots.plot(TT, (+3*V))

Plots.plot!(T[1:570], -S[1:570]*1e3*3389.5)
Plots.plot!(T, +S*1e3*3389.5)

Plots.plot(TT, centerlist[1, :]*1e3*3389.5-Z4[1, 1:10:end])
Plots.plot!(TT, centerlist[1, 1:]*1e3*3389.5-Z4[1, 1:10:end]+S[1:76]*1e3*3389.5)
Plots.plot!(T, centerlist[1, :]*1e3*3389.5-Z4[1, 1:10:5001]-S*1e3*3389.5)
Plots.plot!(T, centerlist[1, :]*1e3*3389.5-Z4[1, 1:10:5001]+S*1e3*3389.5)
Plots.plot!(TT, traj[1, :, 2]-Z4[1, :], linewidth=0.1, legend = false)


S = [inv(Alist[:, :, i])[2, 2] for i=1:1:length(0:1:100)]
Plots.plot(0:1:75, centerlist[2, 1:76]*1e3*3389.5-Z4[2, 1:100:7501])
Plots.plot!(0:1:75, centerlist[2, 1:76]*1e3*3389.5-Z4[2, 1:100:7501]+S[1:76]*1e3*3389.5)
Plots.plot!(0:1:75, centerlist[2, 1:76]*1e3*3389.5-Z4[2, 1:100:7501]-S[1:76]*1e3*3389.5)
Plots.plot!(T, centerlist[2, :]*1e3*3389.5-Z4[1, 1:10:5001]+S*1e3*3389.5)
V = [sqrt(var[2, 2, i]) for i=1:1:length(0:0.01:75)]
Plots.plot!(0:0.01:75, (avg[2, 1:7501]-Z4[2, 1:7501]-3*V))
Plots.plot!(0:0.01:75, (avg[2, 1:7501]-Z4[2, 1:7501]+3*V))

S = [inv(Alist[:, :, i])[7, 7] for i=1:1:length(0:1:100)]
Plots.plot(0:1:75, centerlist[7, 1:76]*1e3*7.0-Z4[8, 1:100:7501])
Plots.plot!(0:1:75, centerlist[7, 1:76]*1e3*7-Z4[8, 1:100:7501]+S[1:76]*1e3*7)
Plots.plot!(0:1:75, centerlist[7, 1:76]*1e3*7-Z4[8, 1:100:7501]-S[1:76]*1e3*7)
V = [sqrt(var[8, 8, i]) for i=1:1:length(0:0.01:75)]
Plots.plot!(0:0.01:75, (avg[8, 1:7501]-Z4[8, 1:7501]-3*V))
Plots.plot!(0:0.01:75, (avg[8, 1:7501]-Z4[8, 1:7501]+3*V))

Plots.plot!(0:0.01:75, traj[8, 1:7501, 43]-Z4[8, 1:7501], linewidth=0.1, legend = false)
