#Entry Simulation

using StaticArrays
using LinearAlgebra
using TrajectoryOptimization

include("vehicle_geometry.jl")
include("aerodynamics.jl")
include("interpolation.jl")
include("dynamics\\entry_vehicle_6dof.jl")

geometry = SphereconeGeometry(70.0*pi/180, 1.325, 0.2, [0.001;0.;-0.189])
coeff_interp = Coeff_interp(tables(geometry,0.0,181.0,1.0,14)...)
