# Models for Propagation

abstract type AbstractModel end

struct ThreeDOF{T} <: AbstractModel
    vehicle::EntryVehicle{T}           # Entry Vehicle (Geometry and Mass Properties)
    environment::Environment{T}        # Environement (Planet and Atmosphere)
    C_D::Vector{T}                     # Drag Coefficient Table
    C_L::Vector{T}                     # Lift Coefficient Table
end

struct SixDOF{T} <: AbstractModel
    vehicle::EntryVehicle{T}           # Entry Vehicle (Geometry and Mass Properties)
    environment::Environment{T}        # Environement (Planet and Atmosphere)
    table_coeff_forces::Array{T, 2}    # Non-dimensionalized Forces Coefficients
    table_coeff_moments::Array{T, 2}   # Non-dimensionalized Moments Coefficients
    table_coeff_damping::Array{T, 2}   # Non-dimensionalised Damping Coefficients
end

function ThreeDOFModelPhoenixMars()
    vehicle = PhoenixVehicle()
    env = MarsEnvironmentJ2()
    table_CD, table_CL = drag_lift_table(vehicle)
    return ThreeDOF{Float64}(vehicle,
                            env,
                            table_CD,
                            table_CL)
end

function SixDOFModelPhoenixMars()
    vehicle = PhoenixVehicle()
    env = MarsEnvironment()
    coeff_force, coeff_moment, coeff_damping = table_aero_spherecone(vehicle)
    return SixDOF{Float64}(vehicle,
                            env,
                            coeff_force,
                            coeff_moment,
                            coeff_damping)
end
