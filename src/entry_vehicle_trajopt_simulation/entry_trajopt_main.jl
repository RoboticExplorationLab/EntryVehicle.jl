#Main file for entry_vehicle with trajectory optimization
using TrajectoryOptimization
using MeshCat
using CoordinateTransformations
using GeometryTypes: GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh
using Colors: RGBA, RGB
using MeshIO
using FileIO
using LinearAlgebra
using ODE
using DifferentialEquations
using Plots
pyplot()

include("aerodynamic_coeff.jl")
include("quaternions.jl")
include("entry_model_trajopt.jl")
include("animate_traj.jl")

####################################
##########Model Definition##########
####################################

T = Float64
n, m = 13, 1

entry_vehicle_model = Model(dyna!, n, m)

#Model definition
model = entry_vehicle_model
n = model.n; m = model.m
model_d = rk4(model)
model_de = discretize_model(model,:DiffEq_ode78) #can take y solver from DiffEq
####################################
##########Initial State#############
####################################

θ = 91*pi/180 #rotation angle about z-axis body
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)

Re = 3389.5
q00 = Q
qff = [1.0; 0.0; 0.0; 0.0]
x0 = zeros(T, n)
x0[1:3] = [3389.5+125; 0.0; 0.0]/Re
x0[4:7] = q00
x0[8:10] = [0.0; 1.0; 0.0]
x0[11:13] = [0.0; 0.0; 0.0]

####################################
##########Target State##############
####################################

M2 = [0.0 0.0 -1.0;
    -sin(θ) cos(θ) 0.0;
     cos(θ) sin(θ) 0.0] #nose down

Q2 = mat2quat(M2)
Q2 = qconj(Q2)

#expected values after 100 seconds and no control inputs
xf = zero(x0)
xf[4:7] = Q2
xf[11:13] = [0.0; 0.0; 0.0]
#here for 550 second
xf[1:3] = [0.9688826924211867;0.08046302598015104; 0.0]
xf[8:10] = [-0.3218613947994733; 0.4350535668863166; 0.0]

#xdot = zeros(n)
#x, u = x0, [0.0]
#evaluate!(xdot,model,x,u)

####################################
##########Cost Function#############
####################################

N = 101
Q = Diagonal(0.1I,n)
R = Diagonal(0.1I,m)
Qf = Diagonal(100.0I,n) #Put float everywhere

obj = LQRObjective(Q, R, Qf, xf, N)

####################################
############Constraints#############
####################################

bnd = BoundConstraint(n, m, u_min = -0.1, u_max = 0.1)
constraints = [bnd] #only bounds constraints on my system, put the []
bnd isa ConstraintSet
constraints isa ConstraintSet
CON = Constraints(constraints, N)

####################################
############Solver Options##########
####################################

max_con_viol = 1.0e-2 #been changed # relax at first and then increase progressively
verbose=true

opts_ilqr = iLQRSolverOptions{T}(verbose=verbose,live_plotting=:off)

opts_al = AugmentedLagrangianSolverOptions{T}(verbose=verbose,
    opts_uncon=opts_ilqr)

# opts_pn = ProjectedNewtonSolverOptions{T}(verbose=verbose,
#     feasibility_tolerance=max_con_viol,
#     solve_type=:feasible)

opts_altro = ALTROSolverOptions{T}(verbose=verbose,
    opts_al=opts_al);

####################################
########Problem Definition##########
####################################

#ADD INITIAL SEQUENCE OF CONTROL

t0 = 0
tf = 200.0
prob = TrajectoryOptimization.Problem(model_d, obj, x0 = x0, xf=xf, constraints = CON, N=N, tf=tf)
prob.dt
plot(prob.U)
rollout!(prob)
plot(prob.X)
#savefig("state")
#prob = Problem(model, obj, x0 = x0, integration=:rk4, N=N, tf=tf)
TrajectoryOptimization.solve!(prob, opts_al)

####################################
########Results Plotting############
####################################

Z = prob.X
U = prob.U

####################################
###########Visualization############
####################################

t_sim = t0:prob.dt:tf
QQ = [0.707107;-0.707107; -0.0; -0.0] #Q_image2model

animate_traj(t_sim, Z)

plot(prob.U)
title!("Control Sequence [N.m]")
xlabel!("time [s]")
#savefig("control_sequence")

plot(prob.X)
