# Predictor-Corrector method (MPC basically) for 3DOF entry vehicle
# See Entry Guidance : A unified method by Ping lu



# 3DOF dynamical model #########################################################

function entry_vehicle_3dof_dynamics_energy!(ẋ::AbstractVector,x::AbstractVector#=,u::AbstractVector=#,params,e::AbstractFloat)
    ## States: X ∈ R^6;
    # r radius (spherical coordinates)
    # θ longitude
    # ϕ latitude
    # γ flight path angle
    # ψ headind angle
    # s range-to-go (radians)

    # Derivatives are with respect to energy

    r, θ, ϕ, γ, ψ, s = x #state vector

    # v magnitude of velocity
    @show(r, e)
    v = sqrt(2*((1/r)-e))

    σ = u[1] #bank angle (control input)

    #Parameters
    env = params[:Env]
    vehicle = params[:Vehicle]
    R_p = env.R_p
    ρ0 = env.ρ0
    H = env.H
    A_ref = params[:A_ref]
    L_ref = params[:L_ref]
    m = vehicle.m

    α = 20 # fixed for now

    #Drag and Lift coeff/forces
    C_D = table_CD[α]
    C_L = table_CL[α]

    h = (r*R_p - R_p) #m
    ρ = exponential_atmosphere(h,ρ0,H)
    V = v*sqrt(9.81*R_p) #real vel
    D = (0.5*ρ*V^2*(A_ref/m)*C_D)/(9.81) #Drag Acceleration
    L = 0.5*ρ*V^2*(A_ref/m)*C_L/9.81

    #1-3 are same as simplified version
    ẋ[1] = sin(γ)/D
    ẋ[2] = (cos(γ)*sin(ψ))/(r*D*cos(ϕ))
    ẋ[3] = (cos(γ)*cos(ϕ))/(r*D)
    ẋ[4] = (1/D)*(L*cos(σ)+((v^2)-(1/r))*cos(γ)/r)
    ẋ[5] = (1/(D*(v)^2))*((L*sin(σ)/(cos(γ)))+((v^2)/r)*cos(γ)*sin(ψ)*tan(ϕ))
    ẋ[6] = -cos(γ)/(r*D)
    return ẋ
end

# energy function ##############################################################
function e(r,v)
    return (1/r)-0.5*v^2
end
