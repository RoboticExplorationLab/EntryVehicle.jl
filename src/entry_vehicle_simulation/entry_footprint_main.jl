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

#Main file for footprint drawing using MeshCat

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model.jl")
include("integration.jl")
include("footprint.jl")
include("animate_footprint.jl")

##########################
###### Integration #######
##########################

#PUT INITIAL POINT IN PARAMETER HERE
Re = 3389.5 #3396.2 [km] Radius of the planet
Δt = 500.0 #seconds - length simulation

function dyn(t, x)
    return dyna(t, x, [0.0]) # torques input should be in km^2*kg*s^(-2) so should be small values
end

pos_end, state_end = footprint()

##########################
###### Visualization #####
##########################

animate_footprint(pos_end)

#=
mean(pos_end[:, 1])
mean(pos_end[:, 2])
mean(pos_end[:, 3])

mean(state_end[:, 8])
mean(state_end[:, 9])
mean(state_end[:, 10]) =#
