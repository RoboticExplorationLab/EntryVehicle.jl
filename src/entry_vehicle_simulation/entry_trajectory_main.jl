using MeshCat
using CoordinateTransformations
using GeometryTypes: GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh
using Colors: RGBA, RGB
using MeshIO
using FileIO
using LinearAlgebra
using ODE
using Plots
using DifferentialEquations
#using ApproxFun
#using Libdl
using PyPlot
#pyplot()
gr()

#This file computes an entry trajectory for the given initial parameters...
#and enables its visualization using MeshCat

####################################
##############Scripts###############
####################################

include("quaternions.jl") #useful functions for rotation
include("aerodynamic_coeff.jl") #compute aerodynamic coefficients
include("entry_model.jl") #entire dynamics model
include("integration.jl") #integration process
include("traj_plots.jl") #plots functions
include("animate_traj.jl") #animation function
include("space_mechanics.jl") #space mechanics coordinates

#Used for MARSGRAM use
include("interpolation.jl")
include("marsgram_wrapper.jl")

########################################
#####Vehicle Geometry - Aerodynamics####
########################################

δ = 70*pi/180
r_cone = 1.3
r_min = 0.2
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

####################################
######Simulation Parameters#########
####################################

Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M) #CHANGE THAT
Q = qconj(Q)
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Δt = 280.0

####################################
######Uncertainty Parameters########
####################################

u = [0.0] #control input - torque about z axis
w = [0.0158*10^9; 0.0; 0.0; 0.0] #first value is the atmospheric density
#vector after that is the wind component

####################################
#####Dynamics - Integration#########
####################################

function dyn(t, x)
    #return dyna_aero(t, x, [0.0], table_CF, table_Cτ) # torques input should be in m^2*kg*s^(-2), use dyna instead for online aerodynamic coefficients
    return dyna_coeffoff(t, x, [0.0])
end

t_sim, Z = integration(dyn, x0, Δt) #dyn or dyna_coeffoff_inplace
t_sim, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)
t_sim4, Z4 = rk4(dyna_coeffoff, x0, [0.0], 0.01, [0.0, Δt]) #really good

####################################
######Visualization - MeshCat#######
####################################

#ROTATION NEEDED FOR MODEL-IMAGE IN VISUALIZATION
#M_model_image = [1.0 0.0 0.0;
#    0.0 0.0 1.0;
#    0.0 -1.0 0.0]
QQ = [0.707107;-0.707107; -0.0; -0.0] #Q_image2model

animate_traj(t_sim, Z)


####################################
#####State Variables - 2D Plots#####
####################################

plot_traj(Z4)
#Plots.savefig("1")
plot_altitude(Z4, t_sim4)
#Plots.savefig("2")
plot_vel(Z, t_sim)
#Plots.savefig("3")
plot_quaternions(Z4, t_sim4)
#savefig("4")
plot_ang_vel(Z4, t_sim4)
#savefig("5")
norm_quaternions(Z4, t_sim4)

plot_traj3D(Z)


####################################
######## Test & Other Plots#########
####################################

Plots.plot(Z[1, :], Z[2, :],legend =false)
ylims!((0.38, 0.405))
xlims!((0.90, 0.95))
title!("Trajectory Sampling - Equatorial Orbits")
xlabel!("X/R")
ylabel!("Y/R")

angles = 0.0:0.01:2*pi
A = [1.0 0.0; 0.0 1.0]
b = [0.0;0.0]
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end

####################################
#######MONTE CARLO SIMULATION#######
####################################

using Distributions
using Random

function generate_samples(x_0, Q, M)
    #M number of samples
    univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    D = Product(univariate_D_vector)
    X_samples = zeros(13, M)
    rand!(D, X_samples)
    return X_samples
end

function prop_MC_entry(X_samples, t_start, t_end, dt)
    n, M = size(X_samples)
    #saveAT = 1.0
    T = t_start:dt:t_end
    traj = zeros(n, length(T), M)
    for i=1:1:M
        @show(i)
        x_ini = X_samples[:, i]
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(dyna_coeffoff, x_ini, [0.0], dt, [t_start, t_end])
        traj[:, :, i] = Z
    end
    return traj
end

x_0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], -3.0, 5.0, 0.0, 0.0, 0.0, 0.0]
Q = Diagonal(0.00001*ones(13))
X_samples = generate_samples(x_0, Q, 100)
t_start = 0.0
dt = 0.01
t_end = 300.0
traj = prop_MC_entry(X_samples, t_start, t_end, dt)

plot_traj2(traj[:, :, 11])
plot_altitude()

function plot_traj_MC(traj)
    n, T, M = size(traj)
    for i=1:1:M
        X = traj[1, :, i]
        Y = traj[2, :, i]
        if i == 1
            Plots.plot(X, Y)
        else
            Plots.plot!(X, Y)
        end
    end
end

plot_traj_MC(traj)


Plots.scatter(traj[1, 1, :], traj[2, 1, :])
