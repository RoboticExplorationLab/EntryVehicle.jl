using MeshCat
using CoordinateTransformations
using GeometryTypes: GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh
using Colors: RGBA, RGB
using MeshIO
using FileIO
using LinearAlgebra
using ODE
using HCubature
using StaticArrays
using WebIO
using Statistics
using DifferentialEquations
using Plots
gr()

#Main file for footprint drawing using MeshCat

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model.jl")
include("integration.jl")
include("footprint.jl")
include("animate_footprint.jl")
include("traj_plots.jl")

##########################
###### Integration #######
##########################

#PUT INITIAL POINT IN PARAMETER HERE
Re = 3389.5 #3396.2 [km] Radius of the planet
Δt = 300.0 #seconds - length simulation

function dyn(t, x)
    return dyna(t, x, [0.0]) # torques input should be in km^2*kg*s^(-2) so should be small values
end

δ = 70*pi/180
r_cone = 1.3
r_min = 0.2
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)


pos_end, state_end = footprint()
trajs_all = footprint_all()

##########################
###### Visualization #####
##########################

animate_footprint(pos_end)

plot_traj(trajs_all[:, :, 1])

#=
mean(pos_end[:, 1])
mean(pos_end[:, 2])
mean(pos_end[:, 3])

mean(state_end[:, 8])
mean(state_end[:, 9])
mean(state_end[:, 10]) =#

##########################
######### Plots ##########
##########################

include("space_mechanics.jl")

function convertt(pos_end)
    pos_geoc = zeros(size(pos_end))
    for i=1:1:36
        pos_ecef = eci2ecef(0.1, pos_end[i, :])
        pos_geoc[i, :] = ecef2geoc(pos_ecef)
        #pos_geod[i, :] = geoc2geod(pos_geoc)
    end
    return pos_geoc
end

pos_geoc = convertt(pos_end)

anim = @animate for j=1:1:36
    Plots.scatter(pos_geoc[1:j, 2], pos_geoc[1:j, 3])
    ylims!((10.0, 11.0))
    xlims!((-0.15, 0.15))
    xlabel!("longitude - degrees")
    ylabel!("latitude - degrees")
    title!("Landing footprint")
end
gif(anim, "footprint.gif", fps = 3)

scatter(pos_geoc[:, 2], pos_geoc[:, 3])
ylims!((10.0, 11.0))
xlabel!("longitude - degrees")
ylabel!("latitude - degrees")
title!("Landing footprint")
savefig("footprint_geoc")
