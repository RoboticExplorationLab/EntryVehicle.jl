################################################################################
############################ DUFFING ###########################################
################################################################################

using LinearAlgebra
using Plots

include("mc.jl")
include("pce.jl")
include("ellipsoid_prop.jl")
include("lin_cov.jl")
include("space_mechanics.jl")
include("traj_plots.jl")

model = duffing
model_lincov = duffing_lincov
p = [-1.0; 1.0; 0.2; 0.1; 1.0]

x0 = [0.1; 0.1]
U0 = [0.01; 0.01]
Q0 = Diagonal(U0.^2) #Q is the covariance matrix
S = [1.0; 1.0]

t0 = 0.0
tf = 25.0
dt_rk = 1e-3
dt_e = 0.1
dt_mc = 0.01
dt_pce = 0.01
dt_lincov = 0.1
scheme = ellipse2points
e_solver = DRN_algo
e_solver = points2ellipse_mosek
k = 2
n = length(x0)
n_scheme = 2*n+1
M_mc = 10000
M_pce = 2000
type = "gaussian"
D = MvNormal
d = 5

#Nominal trajectory
t_sim, Z = rk4(model, x0, p, 0.01, [t0; tf])

#Ellipsoid
Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p, k, e_solver)
#Monte Carlo
T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE(t0, tf, dt_mc, M_pce, type, D, d, x0, Q0, model, p)
#LinCov
T_lin, Q_lin = lincov_prop(x0, Q0, Z, t0, tf, dt_lincov, model_lincov)
################################################################################
############################## Plots ###########################################
################################################################################

#NOMINAL AND CENTERLIST

p02 = Plots.scatter(T_e, centerlist[1, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p01 = Plots.plot!(t_sim, Z[1, :], markersize = 2.0, color = :red, label = "nominal trajectory",  xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
xlabel!("time [s]")
ylabel!("position")
savefig("duff_1.pdf")

p02 = Plots.scatter(T_e, centerlist[2, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p01 = Plots.plot!(t_sim, Z[2, :], markersize = 2.0, color = :red, label = "nominal trajectory",  xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
xlabel!("time [s]")
ylabel!("velocity")
savefig("duff_2.pdf")

p04 = Plots.scatter(centerlist[1, :], centerlist[2, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p03 = Plots.plot!(Z[1, :], Z[2, :], markersize = 2.0, color = :red, label = "nominal trajectory", xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
xlabel!("position")
ylabel!("velocity")
savefig("duff_3.pdf")

#PLots ellipses on the phase portrait
j = 10
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    #B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
    B[:, i] = [cos(angles[i]), sin(angles[i])]
end
ellipse  = (Alist[1:2, 1:2, j]) \ B
#Plots.scatter([centerlist[1, j]], [centerlist[2, j]], markercolor = :black, markersize = 3.0, legend = false)
Plots.plot!(#=centerlist[1, j].+ =#ellipse[1, :]#=centerlist[2, j].+=# ,ellipse[2, :], linewidth = 3.0, linestyle = :dot, label = "t = $j", xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
Plots.scatter!([Z[floor(Int, j), 1]], [Z[floor(Int, j), 2]])
xlabel!("position")
ylabel!("velocity")
savefig("duff_4.pdf")

#OTHER METHODS COMPARISON

vari_e = variance_e(Alist, S)
V_mc = sig(var_MC, S)
V_pce = sig_pce(var_PCE)
V_lincov = sig(Q_lin, S)

j = 2
p1 = Plots.plot(T_MC, 3*V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "MC")
p1 = Plots.plot!(T_MC,  -3*V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "")
p2 = Plots.plot!(T_PCE,-V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p2 = Plots.plot!(T_PCE, V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=600:1:800], linewidth = 0.001, color = :blue, label = "")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=1:1:1], linewidth = 0.001, color = :blue, label = "individual traj",  xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
p4 = Plots.plot!(T_e, -vari_e[j, :], linewidth = 3.0, linestyle = :dash, color = :red, label ="ellipsoids")
p4 = Plots.plot!(T_e, vari_e[j, :], linewidth = 3.0, linestyle = :dash, color = :red, label ="")
p4 = Plots.plot!(T_lin, -V_lincov[j, :], linewidth = 3.0, linestyle = :dash, color = :orange, label ="lincov")
p4 = Plots.plot!(T_lin, V_lincov[j, :], linewidth = 3.0, linestyle = :dash, color = :orange, label ="")
xlabel!("time [s]")
ylabel!("velocity dispersion")
savefig("duff_6.pdf")

################################################################################
############################ VINHS MODEL########################################
################################################################################

model = vinhs_full_model
model_lincov = vinhs_full_model_lincov
p = [45.0*pi/180]

S = [1e3*3389.5; 1.0; 1.0; 1e3*7.00; 1.0; 1.0]
x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
U0 = [50.0; (0.0017); (0.0017); (1.0); (0.0017); (0.0017)]
Q0 = Matrix(Diagonal(U0.^2))

t0 = 0.0
tf = 80.0
dt_rk = 1e-3
dt_e = 1e-1
dt_mc = 0.01
dt_pce = 0.01
dt_lincov = 0.01
scheme = ellipse2points
e_solver = DRN_algo
e_solver = points2ellipse_mosek
n = length(x0)
n_scheme = 2*n+1
M_mc = 10000
M_pce = 2000
type = "gaussian"
D = MvNormal
d = 5

#Nominal trajectory
t_sim, Z = rk4(model, x0, p, 0.01, [t0; tf])

#Ellipsoid
Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p, k, e_solver)
#Monte Carlo
T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE(t0, tf, dt_pce, M_pce, type, D, d, x0, Q0, model, p)
#LinCov
T_lin, Q_lin = lincov_prop(x0, Q0, Z, t0, tf, dt_lincov, model_lincov)
################################################################################
############################## Plots ###########################################
################################################################################


#NOMINAL AND CENTERLIST

center_state_plot(T_e, centerlist)
savefig("vinhs_11.pdf")

p01 = Plots.scatter(t_sim, Z[1, :], markersize = 2.0, color = :blue, label = "nominal trajectory")
p02 = Plots.plot!(T_e, centerlist[1, :], linewidth = 2.0, color = :red, label = "ellipse center")
xlabel!("time [s]")
ylabel!("velocity")
savefig("duffing_2.pdf")

p03 = Plots.scatter(Z[1, :], Z[2, :], markersize = 2.0, color = :blue, label = "nominal trajectory")
p04 = Plots.plot!(centerlist[1, :], centerlist[2, :], linewidth = 2.0, color = :red, label = "ellipse center")
xlabel!("position")
ylabel!("velocity")
savefig("duffing_3.pdf")


#OTHER METHODS COMPARISON

vari_e = variance_e(Alist, S)
V_mc = sig(var_MC)
V_pce = sig_pce(var_PCE)
V_lincov = sig(Q_lin)

j = 6
p1 = Plots.plot!(T_MC, 3*V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "MC")
p1 = Plots.plot!(T_MC,  -3*V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "")
p2 = Plots.plot!(T_PCE,-V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p2 = Plots.plot!(T_PCE, V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "")
p3 = Plots.plot(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=2400:1:2800], linewidth = 0.001, color = :blue, label = "")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=1:1:1], linewidth = 0.001, color = :blue, label = "individual traj",  xtickfont = font(9), xguidefontsize=15, ytickfont = font(7), yguidefontsize=13)
p4 = Plots.plot!(T_e, -vari_e[j, :], linewidth = 4.0, linestyle = :dash, color = :red, label ="ellipsoids")
p4 = Plots.plot!(T_e, vari_e[j, :], linewidth = 4.0, linestyle = :dash, color = :red, label ="")
p4 = Plots.plot!(T_lin, -V_lincov[j, :], linewidth = 3.0, linestyle = :dash, color = :orange, label ="lincov")
p4 = Plots.plot!(T_lin, V_lincov[j, :], linewidth = 3.0, linestyle = :dash, color = :orange, label ="", legend = false )
ylabel!("heading dispersion [rad]")
xlabel!("time [s]")
savefig("3dof_66.pdf")

uncertainty(Alist)

rank(XX[:, :, 6000])
cond(XX[:, :, 2], 2)

cond(Alist[:, :, 1], 2)

X_lims(XX[:, :, 6000])


################################################################################
############################ 6 DOF MODEL########################################
################################################################################

#Ellipses
v_eci = [-1.6; 6.8; 0.0001]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β);-sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q)
#v_eci = [-0.001; 1.0; 0.0001]*1e3
x0 = [(3389.5+125)*1e3; 0.0; 50.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]
x0_12_e = [(3389.5+125)*1e3; 0.0; 50.0; 0.00; 0.00; 0.00; v_eci; 0.0; 0.0; 0.0]#because the center with respect to himself is 0
S = [3389.5*1e3;3389.5*1e3;3389.5*1e3; 1.0;1.0;1.0;7.00*1e3;7.00*1e3;7.00*1e3;1.0;1.0;1.0]
U0 = [50.0;50.0;1.0;0.005;0.005;0.005;1e-3;1e-3;1e-3;1e-6;1e-6;1e-6]
A0, b0 = ini(x0_12_e, U0, S)

cond(A0, 2) #really bad need another way of sclaing things
A0

#MC
e = quat2euler(Q)
x0_12 = [(3389.5+125)*1e3; 0.0; 50.0; e; v_eci; 0.0; 0.0; 0.0]
Q0_12 = Matrix(Diagonal(U0.^2))

model = dyna_coeffoff_COM_on_axis
model_lincov =  dyna_coeffoff_COM_on_axis_lincov
t0 = 0.0
tf = 60.0
dt_rk = 1e-2
dt_e = 1.0
dt_mc = 0.01
dt_pce = 0.1
dt_lincov = 0.01
n = length(x0)
M_mc = 1000
M_pce = 20
type = "gaussian"
D = MvNormal
d = 6

p = [0.0]

u0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
R = entry3DOF2ecef(u0)
x00 = [R[1:3]; Q[1]; Q[2]; Q[3]; Q[4]; R[4:6]; 0.0; 0.0; 0.0]
x00_12 = [R[1:3]; zeros(3); R[4:6]; 0.0; 0.0; 0.0]
A0, b0 = ini(x00_12, U0, S)

#Nominal trajectory
t_sim, Z = rk4(model, x0, p, 0.01, [t0; 180])
Z_12 = point13_12_mc_full(Z)

plot_total_attack_angle(Z, t_sim)
savefig("AoA.pdf")

plot_altitude(Z, t_sim)
savefig("altitude.pdf")

plot_ang_vel(Z, t_sim)
savefig("ang_vel.pdf")

plot_entry_profile(Z, t_sim)
savefig("entry_profile.pdf")

plot_mach_number_altitude(Z, t_sim)
savefig("mach_alt.pdf")


Plots.plot(Z_13[13, :])

#Ellipsoid
Alist, blist, centerlist, XX, ref ,T_e = ellipse_propagation(A0, b0, t0, tf, dt_e, scale_down_13(Z[:, 1], S), model, p)
#Monte Carlo

a = 1.0

T_MC, traj_MC, m_MC, var_MC = simu_mc_6dof(x0_12, Q0_12, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE_6dof(t0, tf, dt_pce, M_pce, type, D, d, x0_12, Q0_12, model, p)
#LinCov
T_lin, Q_lin = lincov_prop(x0_12, Q0_12, Z_12, t0, tf, dt_lincov, model_lincov)

#NOMINAL AND CENTER TRAJECTORIES
center = transform_centerlist(centerlist, ref, S) #now in dim 13 with the right scaling

l = 13
Plots.scatter(T_e, center[l, :])
Plots.plot!(t_sim, Z[l,:])

p1 = Plots.scatter(T_e, centerlist[1, :])
p2 = Plots.scatter(T_e, centerlist[2, :])
p3 = Plots.scatter(T_e, centerlist[3, :])
p4 = Plots.scatter(T_e, centerlist[4, :])
p5 = Plots.scatter(T_e, centerlist[5, :])
p6 = Plots.scatter(T_e, centerlist[6, :])
p7 = Plots.scatter(T_e, centerlist[7, :])
p8 = Plots.scatter(T_e, centerlist[8, :])
p9 = Plots.scatter(T_e, centerlist[9, :])
p10 = Plots.scatter(T_e, centerlist[10, :])
p11 = Plots.scatter(T_e, centerlist[11, :])
p12 = Plots.scatter(T_e, centerlist[12, :])

Plots.plot(p1, p2, p3, layout = (1, 3))


#COMPARISON UNCERTAINTY
vari_e = variance_e(Alist, S)
V_mc = sig(var_MC)
V_pce = sig_pce(var_PCE)
V_lincov = sig(Q_lin)

j = 4
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :]+ V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "MC")
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :]+ V_mc[j, :], linewidth = 2.0)
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :]- V_mc[j, :], linewidth = 2.0, label = "")
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :]+ V_mc[j, :], linewidth = 2.0, legend = :topleft)
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :]- V_mc[j, :], linewidth = 2.0, label = "")
p1 = Plots.plot!(T_MC, m_MC[j, :] -Z_12[j, :] -V_mc[j, :], linestyle = :dot, linewidth = 3.0, color = :black, label = "")
p2 = Plots.plot!(T_PCE,-V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p2 = Plots.plot!(T_PCE, V_pce[j, :], linestyle = :dash, linewidth = 2.0, color = :green, label = "")
p3 = Plots.plot(T_MC, [traj_MC[j, :, i]-Z_12[j, :] for i=600:1:800], linewidth = 0.001, color = :blue, label = "")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z_12[j, :] for i=1:1:1], linewidth = 0.001, color = :blue, label = "individual traj",  xtickfont = font(9), xguidefontsize=15, ytickfont = font(8), yguidefontsize=13)
p4 = Plots.plot!(T_e, -vari_e[j, :], linewidth = 3.0, color = :red, label ="ellipsoids")
p4 = Plots.plot!(T_e, +vari_e[j, :], linewidth = 3.0, color = :red, label ="")
p4 = Plots.plot!(T_lin, -V_lincov[j, :], linewidth = 3.0, color = :orange, label ="lincov")
p4 = Plots.plot!(T_lin, V_lincov[j, :], linewidth = 3.0, color = :orange, label ="", legend = :topleft)
ylabel!("e x [rad]")
xlabel!("time [s]")
savefig("6dof_444.pdf")
ylims!((-0.03,0.03))
ylabel!("psi dispersion [rad]")
xlabel!("time [s]")
savefig("vinhs_compare_6.pdf")






X = [1.035 1.03497 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498; 0.00802473 0.00802473 0.00803948 0.00800998 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473; 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48991e-5 1.48401e-5 1.48697e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48697e-5 1.48695e-5 1.48696e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554; 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413; 1.43241e-5 1.43248e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43245e-5 1.43567e-5 1.42922e-5 1.43238e-5 1.43244e-5 1.43242e-5 1.43246e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43387e-5 1.43101e-5 1.4339e-5 1.43098e-5 1.43242e-5 1.43245e-5 1.43244e-5 1.43244e-5 1.43244e-5; -27915.1 -28203.9 -28058.6 -28059.7 -28059.2 -28059.2 -54758.0 -1056.1 -31359.8 -30101.5 -27936.8 -28181.2 -28059.2 -28059.1 -28059.2 -28059.1 -28059.2 -28059.1 59795.6 1.15791e5 -29027.7 -28686.3 -28060.1 -28058.1 -28059.2; -1.66149 -1.67918 -1.67173 -1.6689 -1.67031 -1.67031 -2.05381 -0.17184 -17.5735 13.3792 -1.59774 -1.7436 -1.67072 -1.66991 -1.67042 -1.67021 -1.67032 -1.67031 -1.82819 -1.44326 91.757 -95.1843 -1.62103 -1.7196 -1.67031; -5.87788 -5.94103 -5.90923 -5.90951 -5.90937 -5.90937 -5.97552 -5.83186 -5.92979 -5.92207 -5.90903 -5.9097 -5.90937 -5.90937 -5.90938 -5.90937 -5.90937 -5.90937 -5.94055 -5.87516 -5.91493 -5.91297 1.0e16 -1.0e16 -5.90937]

rank(vcat(X, ones(25)'))


X = [1.03528 1.03469 1.03438 1.03458 1.03498 1.03498 1.03468 1.03498 1.03498 1.03498 1.03498 1.03498 1.03616 1.0338 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498; 0.00802474 0.00802473 0.00831976 0.00772971 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802474 0.00802473 0.00920484 0.00684463 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473; 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 0.000309896 -0.000280157 1.48697e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48695e-5 0.00119498 -0.00116524 1.48696e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5; -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.0876957 -0.373413 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554; 0.971414 0.971412 0.971413 0.971414 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971415 0.97141 1.11427 0.828559 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413; 1.43177e-5 1.43318e-5 1.43244e-5 1.43244e-5 1.37591e-5 1.48897e-5 1.43567e-5 1.42922e-5 1.43238e-5 1.43244e-5 1.43242e-5 1.43246e-5 1.43064e-5 1.43465e-5 1.43419e-5 1.43056e-5 0.14287 -0.142842 1.4339e-5 1.43098e-5 1.43242e-5 1.43245e-5 1.43244e-5 1.43244e-5 1.43244e-5; -25309.8 -31093.5 -28048.2 -28069.4 -28056.9 -28060.6 -54758.0 -1056.1 -31359.8 -30101.5 -27936.8 -28181.2 -40262.5 -54551.8 -59517.0 -32429.2 -7.44513e5 6.51551e5 59795.6 -1.15791e5 -29027.7 -28686.3 -28060.1 -28058.1 -28059.2; -1.50214 -1.85649 -1.69859 -1.64198 -1.67027 -1.67031 -2.05381 -0.17184 -17.5735 13.3792 -1.59774 -1.7436 -526.069 883.295 -182.306 138.228 4.34169 -9.77991 -1.82819 -1.44326 91.757 -95.1843 -1.62103 -1.7196 -1.67031; -5.31041 -6.57585 -5.90654 -5.91203 -5.90927 -5.90929 -5.97552 -5.83186 -5.92977 -5.92208 -5.90903 -5.90971 -4.20089 -7.93326 -8.01631 -4.53088 -5.62155 -5.79802 -5.94058 -5.87516 -5.91493 -5.91312 1.0e16 -1.0e16 -5.90937]


X = [1.035 1.03497 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498 1.03498; 0.00802473 0.00802473 0.00803948 0.00800998 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473 0.00802473; 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 2.96209e-5 1.1826e-7 1.48697e-5 1.48695e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48708e-5 1.48684e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5 1.48696e-5; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230555 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554 -0.230554; 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413 0.971413; 1.43241e-5 1.43248e-5 1.43244e-5 1.43244e-5 1.42961e-5 1.43527e-5 1.43567e-5 1.42922e-5 1.43238e-5 1.43244e-5 1.43242e-5 1.43246e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.44673e-5 1.41816e-5 1.43259e-5 1.4323e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43244e-5 1.43244e-5; -27915.1 -28203.9 -28058.6 -28059.7 -28059.1 -28059.2 -54758.0 -1056.1 -31359.8 -30101.5 -27936.8 -28181.2 -28059.2 -28059.1 -28059.2 -28059.1 -28060.0 -28058.4 -19279.2 -36837.9 -28085.6 -28049.9 -28059.3 -28059.1 -28059.2; -1.66149 -1.67918 -1.67173 -1.6689 -1.67031 -1.67032 -2.05381 -0.17184 -17.5735 13.3792 -1.59774 -1.7436 -1.67072 -1.66991 -1.67042 -1.67021 -1.67033 -1.6703 -1.68859 -1.65136 7.70136 -11.0429 -1.66539 -1.67524 -1.67031; -5.87788 -5.94103 -5.90923 -5.90951 -5.90937 -5.90937 -5.97552 -5.83186 -5.92979 -5.92207 -5.90903 -5.9097 5.90937 -5.90937 -5.90938 -5.90937 -5.90937 -5.90937 -5.91262 -5.90609 -5.90952 -5.90929 1.0e15 -1.0e15 -5.90937]

X[:, 3] == X[:, 10]

rank(X, atol = 1e-40, rtol = 1e-40)

svd(X).S

T = BigFloat
A = zeros(T, 3)

S = [8.70472e-10 4.59238e-22 1.21203e-14 0.0 0.0 0.0 -5.32907e-15 7.10543e-15 3.25261e-19 -2.84425e-20 -2.84425e-20 -2.84425e-20; -1.35412e-22 8.7042e-10 2.64518e-20 0.0 0.0 0.0 1.72686e-22 9.74232e-22 1.58653e-26 6.57614e-35 6.57614e-35 6.57614e-35; 1.21201e-14 2.64518e-20 8.7042e-10 0.0 0.0 0.0 -5.42101e-20 0.0 3.30872e-24 -2.15702e-25 -2.15702e-25 -2.15702e-25; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 3.55271e-15 8.33429e-22 -5.42101e-20 0.0 0.0 0.0 7.99361e-14 0.0 2.71051e-20 -4.73146e-22 -4.73146e-22 -4.73146e-22; 7.10543e-15 1.35525e-22 0.0 0.0 0.0 0.0 0.0 8.52651e-14 -2.1684e-19 -1.0842e-21 -1.0842e-21 -1.0842e-21; 0.0 2.35778e-26 -3.30872e-24 0.0 0.0 0.0 2.71051e-20 -4.33681e-19 8.16327e-14 2.77033e-25 2.77033e-25 2.77033e-25; -1.32338e-20 -2.37817e-35 3.00089e-26 0.0 0.0 0.0 -2.26536e-21 -1.17825e-20 1.29102e-25 4.0e-8 0.0 0.0; -1.32338e-20 -2.37817e-35 3.00089e-26 0.0 0.0 0.0 -2.26536e-21 -1.17825e-20 1.29102e-25 0.0 4.0e-8 0.0; -1.32338e-20 -2.37817e-35 3.00089e-26 0.0 0.0 0.0 -2.26536e-21 -1.17825e-20 1.29102e-25 0.0 0.0 4.0e-8]
