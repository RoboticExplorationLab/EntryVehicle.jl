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
pyplot()

#This file computes an entry trajectory for the given initial parameters...
#and enables its visualization using MeshCat

####################################
##############Scripts###############
####################################

include("quaternions.jl") #useful functions for rotation
include("aerodynamic_coeff.jl") #compute aerodynamic coefficients
include("entry_model_dis.jl") #entire dynamics model
include("propagation.jl")
include("traj_plots.jl") #plots functions
include("animate_traj.jl") #animation function


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
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 300, 0.0]
#x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 2.0, 0.0, 0.0, 0.0, 0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]
Δt = 500 #length simulation
dt = 0.001

####################################
#####Dynamics - Integration#########
####################################

t_sim, Z = propagation(x0, Δt, dt)

####################################
######Visualization - MeshCat#######
####################################

#ROTATION NEEDED FOR MODEL-IMAGE IN VISUALIZATION
#M_model_image = [1.0 0.0 0.0;
#    0.0 0.0 1.0;
#    0.0 -1.0 0.0]
QQ = [0.707107;-0.707107; -0.0; -0.0] #Q_image2model

#Animate Trajectory
animate_traj(t_sim[1:1000:end], Z[:, 1:1000:end])

####################################
#####State Variables - 2D Plots#####
####################################

plot_traj(Z)
#savefig("1")
plot_altitude(Z, t_sim)
#savefig("2")
plot_vel(Z, t_sim)
#savefig("3")
plot_quaternions(Z)
#savefig("4")
plot_ang_vel(Z)
#savefig("5")
norm_quaternions(Z, t_sim)

plot_temp(Z, t_sim)
savefig("Temperature_1")

plot_acc(Z, t_sim)
