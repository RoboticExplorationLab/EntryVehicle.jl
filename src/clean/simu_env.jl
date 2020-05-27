#Entry Simulation
using LinearAlgebra

# TO DO:
# Inteface with the right atmosphere and gravitation functions properly
# Use Interpolations functions for both models
# Check All Aerodynamics Coefficients
# Add proper Reference Frames
# Add proper Time Reference Frame (with rotating planet)
# CHeck Units Agreement on Model Parameters/Dyamics/Coefficients

include("Vehicles.jl")
include("Aerodynamics.jl")
include("Atmospheres.jl")
include("Environments.jl")
include("Models.jl")
include("Integration.jl")
include("dynamics\\Quaternions.jl")
include("dynamics\\entry_vehicle_3dof.jl")
include("SpaceMechanics.jl")
#include("Interpolation.jl")

# Define a Vehicle
#phoenix_vehicle = PhoenixVehicle()
# Define a Planetary Environment (contains atmosphere and )
#environment = MarsEnvironment()

# Define a Dynamics Model (3 DOF, 6 DOF)
model = ThreeDOFModelPhoenixMars()
dynamics = entry_vehicle_3dof_dynamics

# Define Initial Conditions
v_eci = [-1.6; 6.8; 0.0001]*1e3
r_eci = [(3389.5+125)*1e3; 0.0; 50.0]
r_entry3DOF = ecef2entry3DOF([r_eci; v_eci])
γ = r_entry3DOF[5]*180/pi    # -13 degrees in flight path angle entry
y_0 = r_entry3DOF
dt = 1e-2
t_span = [0.0, 100.0]

# Integrate Model
T_sim, Y_sim = rk4(dynamics, y_0, model, dt, t_span)
Plots.plot(T_sim, Y_sim[1, :])  # Plot Altitude Decrease
