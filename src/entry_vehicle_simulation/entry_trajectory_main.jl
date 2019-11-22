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
x0 = [(3389.5+900)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], -3.0, 5.0, 0.0, 0.0, 0.0, 0.0]
Δt = 1200.0

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

plot_traj(Z)
#Plots.savefig("1")
plot_altitude(Z, t_sim)
#Plots.savefig("2")
plot_vel(Z, t_sim)
#Plots.savefig("3")
plot_quaternions(Z, t_sim)
#savefig("4")
plot_ang_vel(Z, t_sim)
#savefig("5")
norm_quaternions(Z, t_sim)

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
