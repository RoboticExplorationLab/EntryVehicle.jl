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
pyplot()
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


####################################
######Simulation Parameters#########
####################################

a = 1

δ = 70*pi/180
r_cone = 1.3
r_min = 0.2
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_cone, r_min, r_G)

a = 1

function simulate(θ)

    Re = 3389.5
    #θ = 91*pi/180 #rotation angle about z-axis
    M = [-sin(θ) cos(θ) 0.0;
         0.0 0.0 1.0;
         cos(θ) sin(θ) 0.0]
    Q = mat2quat(M) #CHANGE THAT
    Q = qconj(Q)
    #x0 = [(3389.5+125)/Re, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    Δt = 280 #length simulation

####################################
#####Dynamics - Integration#########
####################################

    #new part for offline aerodynamic coefficients computation
    δ = 70*pi/180
    r_cone = 1.3
    r_min = 0.2
    r_G = [0.2; 0.0; 0.3]
    table_CF, table_Cτ = table_aero_spherecone(δ, r_cone, r_min, r_G)

    function dyn(t, x)
        #return dyna_aero(t, x, [0.0], table_CF, table_Cτ) # torques input should be in m^2*kg*s^(-2), use dyna instead for online aerodynamic coefficients
        return dyna_coeffoff(t, x, [0.0])
    end

    t_sim, Z = integration2(dyna_coeffoff_inplace!, x0, Δt) #dyn or dyna_coeffoff_inplace

    return t_sim, Z
end

#list theta

function multiple()
    Θ = 91:1:100
    Θ = Θ*pi/180
    ZZ = zeros(13, 2801, length(Θ))
    for i=1:length(Θ)
        t_sim, Z = simulate(Θ[i])
        ZZ[:, :, i] = Z
    end
    return t_sim, ZZ
end

function plot_multiple()
    Θ = 91:1:100
    for i=1:length(Θ)
        plot_traj2(ZZ[:, :, i])
    end
end



t_sim, ZZ = multiple()
plot_multiple()
plot_traj(ZZ[:, :, 1])
plot_traj2(ZZ[:, :, 2])
plot_traj2(ZZ[:, :, 10])

plot_altitude(ZZ[:, :, 1], t_sim)
plot_altitude(ZZ[:, :, 2], t_sim)
plot_altitude(ZZ[:, :, 10], t_sim)

plot_quaternions(ZZ[:, :, 1])
plot_quaternions(ZZ[:, :, 2])
plot_quaternions(ZZ[:, :, 10])

####################################
######Visualization - MeshCat#######
####################################

#ROTATION NEEDED FOR MODEL-IMAGE IN VISUALIZATION
#M_model_image = [1.0 0.0 0.0;
#    0.0 0.0 1.0;
#    0.0 -1.0 0.0]
QQ = [0.707107;-0.707107; -0.0; -0.0] #Q_image2model

#Animate Trajectory
animate_traj(t_sim, ZZ[:, :, 10])

#MeshCat.convert_frames_to_video("C:\\Users\\33645\\Downloads\\meshcat_1565028345425.tar", overwrite=false) #not working so far
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
