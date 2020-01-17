################################################################################
############################ DUFFING ###########################################
################################################################################

include("mc.jl")
include("pce.jl")
include("ellipsoid_prop.jl")
include("space_mechanics.jl")
include("traj_plots.jl")

model = duffing
p = [-1.0; 1.0; 0.2; 0.1; 1.0]

x0 = [0.1; 0.1]
U0 = [0.01; 0.01]
Q0 = Diagonal(U0.^2)
S = [1.0; 1.0]

t0 = 0.0
tf = 25.0
dt_rk = 0.01
dt_e = 0.1
dt_mc = 0.01
dt_pce = 0.01
scheme = ellipse2points4
n = length(x0)
n_scheme = 2*n
M_mc = 10000
M_pce = 2000
type = "gaussian"
D = MvNormal

#Nominal trajectory
t_sim, Z = rk4(model, x0, p, dt_rk, [t0; tf])

#Ellipsoid
Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p)
#Monte Carlo
T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE(t0, tf, dt, M_pce, type, D, d, x0, Q0/9, model, p)
#LinCov

################################################################################
############################## Plots ###########################################
################################################################################

#NOMINAL AND CENTERLIST

p02 = Plots.scatter(T_e, centerlist[1, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p01 = Plots.plot!(t_sim, Z[1, :], markersize = 2.0, color = :red, label = "nominal trajectory")
xlabel!("time [s]")
ylabel!("position")
savefig("duffing_1.pdf")

p02 = Plots.scatter(T_e, centerlist[2, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p01 = Plots.plot!(t_sim, Z[2, :], markersize = 2.0, color = :red, label = "nominal trajectory")
xlabel!("time [s]")
ylabel!("velocity")
savefig("duffing_2.pdf")

p04 = Plots.scatter(centerlist[1, :], centerlist[2, :], linewidth = 2.0, color = :blue, label = "ellipse center")
p03 = Plots.plot!(Z[1, :], Z[2, :], markersize = 2.0, color = :red, label = "nominal trajectory")
xlabel!("position")
ylabel!("velocity")
savefig("duffing_3.pdf")

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
Plots.plot!(#=centerlist[1, j].+ =#ellipse[1, :]#=centerlist[2, j].+=# ,ellipse[2, :], linewidth = 2.0, linestyle = :dot, label = "$j")
Plots.scatter!([Z[floor(Int, j), 1]], [Z[floor(Int, j), 2]])
xlabel!("position")
ylabel!("velocity")
savefig("duffing_4.pdf")

#OTHER METHODS COMPARISON

vari_e = variance_e(Alist)
V_mc = sig(var_MC)
V_pce = sig_pce(var_PCE)

j = 2
p1 = Plots.plot(T_MC, [3*V_mc[j, :], -3*V_mc[j, :]], linestyle = :dot, linewidth = 3.0, color = :black, label = "MC", legend = false)
p2 = Plots.plot!(T_PCE, [3*V_pce[j, :], -3*V_pce[j, :]], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=400:1:500], linewidth = 0.001, color = :blue, label = "individual traj")
p4 = Plots.plot!(T_e, [vari_e[j, :], -vari_e[j, :]], linewidth = 3.0, linestyle = :dash, color = :red, label ="ellipsoids")
xlabel!("time [s]")
ylabel!("velocity dispersion")
savefig("duffing_6.pdf")

################################################################################
############################ VINHS MODEL########################################
################################################################################

model = vinhs_full_model
p = [45.0*pi/180]

S = [3389.53*1e3; 1.0; 1.0; 1e3*7.00; 1.0; 1.0]
x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
U0 = [50.0; (0.0017); (0.0017); (1.0); (0.0001); (0.0001)]
Q0 = Matrix(Diagonal(U0.^2))

t0 = 0.0
tf = 80.0
dt_rk = 0.01
dt_e = 1.0
dt_mc = 0.01
dt_pce = 0.01
scheme = ellipse2points
n = length(x0)
n_scheme = 2*n+1
M_mc = 10000
M_pce = 2000
type = "gaussian"
D = MvNormal
d = 5

#Nominal trajectory
t_sim, Z = rk4(model, x0, p, dt_rk, [t0; tf])

#Ellipsoid
Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p)
#Monte Carlo
T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE(t0, tf, dt_pce, M_pce, type, D, d, x0, Q0, model, p)
#LinCov

################################################################################
############################## Plots ###########################################
################################################################################


#NOMINAL AND CENTERLIST

center_state_plot(T_e, centerlist)
savefig("vinhs_11.pdf")

p01 = Plots.scatter(t_sim, Z[2, :], markersize = 2.0, color = :blue, label = "nominal trajectory")
p02 = Plots.plot!(T_e, centerlist[2, :], linewidth = 2.0, color = :red, label = "ellipse center")
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

j = 6
p1 = Plots.plot(T_MC, [3*V_mc[j, :], -3*V_mc[j, :]], linestyle = :dot, linewidth = 2.0, color = :black, label = "MC", legend = false)
p2 = Plots.plot!(T_PCE, [V_pce[j, :], -V_pce[j, :]], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z[j, :] for i=900:1:1000], linewidth = 0.001, color = :blue, label = "individual traj")
p4 = Plots.plot!(T_e, [vari_e[j, :], -vari_e[j, :]], linewidth = 3.0, linestyle = :dash, color = :red, label ="ellipsoids")
ylabel!("psi dispersion [rad]")
xlabel!("time [s]")
savefig("vinhs_compare_6.pdf")




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
x0_12_e = [(3389.5+125)*1e3; 0.0; 50.0; 0.0; 0.0; 0.0; v_eci; 0.0; 0.0; 0.0]#because the center with respect to himself is 0
S = [(1e3*3389.5);(1e3*3389.5); (1e3*3389.5); 1.0; 1.0 ;1.0;(1e3*7.00);(1e3*7.00);(1e3*7.00);1.0;1.0;1.0]
U0 = [50.0;50.0;1.0;0.005;0.005;0.005;1e-3;1e-3;1e-4;1e-40;1e-40;1e-40]
A0, b0 = ini(x0_12_e, U0, S)


#MC
e = quat2euler(Q)
x0_12 = [(3389.5+125)*1e3; 0.0; 50.0; e; v_eci; 0.0; 0.0; 0.0]
Q0_12 = Matrix(Diagonal(U0.^2))


model = dyna_coeffoff_COM_on_axis
t0 = 0.0
tf = 100.0
dt_rk = 0.01
dt_e = 1.0
dt_mc = 0.01
dt_pce = 0.01
n = length(x0)
M_mc = 2000
M_pce = 2000
type = "gaussian"
D = MvNormal
d = 5

p = [0.0]

u0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
R = entry3DOF2ecef(u0)
x00 = [R[1:3]; Q[1]; Q[2]; Q[3]; Q[4]; R[4:6]; 0.0; 0.0; 0.0]
x00_12 = [R[1:3]; zeros(3); R[4:6]; 0.0; 0.0; 0.0]
A0, b0 = ini(x00_12, U0, S)

#Nominal trajectory
t_sim, Z_13 = rk4(model, x00, p, dt_rk, [t0; tf])
Z_12 = point13_12_mc_full(Z_13)

plot_total_attack_angle(Z_13, t_sim)
Plots.plot(Z_13[13, :])


#Ellipsoid
Alist, blist, centerlist, XX, ref ,T_e = ellipse_propagation(A0, b0, t0, tf, dt_e, scale_down_13(x0), model, p)
#Monte Carlo
T_MC, traj_MC, m_MC, var_MC = simu_mc_6dof(x0_12, Q0_12, M_mc, model, p, t0, tf, dt_mc)
#PCE
T_PCE, m_PCE, var_PCE = simulation_PCE_6dof(t0, tf, dt_pce, M_pce, type, D, d, x0_12, Q0_12, model, p)
#LinCov


#NOMINAL AND CENTER TRAJECTORIES
center = transform_centerlist(centerlist, ref, S)


#COMPARISON UNCERTAINTY
vari_e = variance_e(Alist, S)
V_mc = sig(var_MC)
V_pce = sig_pce(var_PCE)

j = 3
p1 = Plots.plot(T_MC, [3*V_mc[j, :], -3*V_mc[j, :]], linestyle = :dot, linewidth = 2.0, color = :black, label = "MC", legend = false)
p2 = Plots.plot!(T_PCE, [V_pce[j, :], -V_pce[j, :]], linestyle = :dash, linewidth = 2.0, color = :green, label = "PCE")
p3 = Plots.plot!(T_MC, [traj_MC[j, :, i]-Z_13[j+1, :] for i=1:1:200], linewidth = 0.001, color = :blue, label = "individual traj")
p4 = Plots.plot(T_e, [vari_e[j, :], -vari_e[j, :]], linewidth = 3.0, linestyle = :dash, color = :red, label ="ellipsoids")
ylims!((-1e3,1e3))
ylabel!("psi dispersion [rad]")
xlabel!("time [s]")
savefig("vinhs_compare_6.pdf")
