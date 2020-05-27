#3 dof entry vehicle dynamics based on Vinh's equations

#TO DO
#add break point when reaching the ground
#see how to deal with control input params or u ? and with respect to time.
#need time in the function
#Interpolate CD and CL

# In place function
function entry_vehicle_3dof_dynamics!(ẋ::AbstractVector,x::AbstractVector,u::AbstractVector,params)
    ## States: X ∈ R^6;
    # r radius (spherical coordinates)
    # θ longitude
    # ϕ latitude
    # v velocity magnitude
    # γ flight path angle
    # ψ headind angle

    r, θ, ϕ, v, γ, ψ = x #state vector

    σ = u[1] #bank angle (control input)

    #Parameters
    req = params[:R_p]
    Re = req
    J2 = params[:J2]
    ω = params[:ω_p][end]
    μ = params[:μ_p]
    m = params[:m]

    #Gravity terms
    g_r = -μ/(r^2)*(1-1.5*J2*((req/r)^2)*(3*(sin(ϕ)^2)-1))
    g_ϕ = -3*J2*(μ/(r^2))*((req/r)^2)*sin(ϕ)*cos(ϕ)

    #Coriolis terms
    C_γ = 2*ω*cos(ψ)*cos(ϕ) #Coriolis
    C_ψ = 2*ω*(tan(γ)*sin(ψ)*cos(ϕ)-sin(ϕ)) #Coriolis

    #Transport terms
    Γ_v = r*(ω^2)*cos(ϕ)*(sin(γ)*cos(ϕ)-cos(γ)*sin(ψ)*sin(ϕ))
    Γ_γ = (r/v)*(ω^2)*cos(ϕ)*(sin(γ)*sin(ψ)*sin(ϕ)+cos(γ)*cos(ϕ))
    Γ_ψ = -(r/v)*(ω^2)*(cos(ψ)*sin(ϕ)*cos(ϕ))/(cos(γ))

    #Attack angle profile
    if  t <= 50.0
        α = 50.0
    else
        α = -t+110
    end

    #α = exp(-0.1*t)*cos(t) +15
    α = floor(Int, α)

    #Drag and Lift coeff/forces
    C_D = table_CD[α]
    C_L = table_CL[α]

    h = r - Re #m
    D = SVector(0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_D) #Drag Acceleration
    L = SVector(0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_L)

    #1-3 are same as simplified version
    ẋ[1] = v*sin(γ)
    ẋ[2] = (v/r)*(cos(γ)*cos(ψ)/cos(ϕ))
    ẋ[3] = (v/r)*cos(γ)*sin(ψ)
    ẋ[4] = -D+g_r*sin(γ)+g_ϕ*cos(γ)*sin(ψ) + Γ_v
    ẋ[5] = (1/v)*(L*cos(σ)+g_r*cos(γ)-g_ϕ*sin(γ)*sin(ψ))+(v/r)*cos(γ) + C_γ + Γ_γ
    ẋ[6] = -(1/(v*cos(γ)))*(L*sin(σ)-g_ϕ*cos(ψ)) -(v/r)*cos(γ)*cos(ψ)*tan(ϕ) +C_ψ + Γ_ψ

    return ẋ
end

# Not in place function
function entry_vehicle_3dof_dynamics(t::T, x::AbstractVector{T}, u::AbstractVector{T}, model::AbstractModel) where {T}
    ## States: X ∈ R^6;
    # r radius (spherical coordinates)
    # θ longitude
    # ϕ latitude
    # v velocity magnitude
    # γ flight path angle
    # ψ headind angle

    # Model contains everything for:
          # Vehicle (Geometry and Mars Properties)
          # Environment (Orbiting Body, GravityField, Atmosphere)
          # Aerodynamics (depending on dynamic model)

    # Initialize output
    ẋ = zeros(length(x))

    # Unpack state's components
    r, θ, ϕ, v, γ, ψ = x #state vector

    # Unpack control's components
    σ = u[1] #bank angle (control input)

    # Unpack Model Parameters
    environment  = model.environment
    vehicle = model.vehicle
    req = environment.R
    Re = req
    J2 = environment.gravity.J2
    ω = environment.Ω
    μ = environment.μ
    m = vehicle.m
    table_CD = model.C_D
    table_CL = model.C_L
    ρ0 = environment.atmosphere.ρ0
    H = environment.atmosphere.H
    A_ref = get_A_ref(vehicle.geometry)
    L_ref = get_L_ref(vehicle.geometry)

    #Gravity terms
    g_r = -μ/(r^2)*(1-1.5*J2*((req/r)^2)*(3*(sin(ϕ)^2)-1))
    g_ϕ = -3*J2*(μ/(r^2))*((req/r)^2)*sin(ϕ)*cos(ϕ)

    #Coriolis terms
    C_γ = 2*ω*cos(ψ)*cos(ϕ) #Coriolis
    C_ψ = 2*ω*(tan(γ)*sin(ψ)*cos(ϕ)-sin(ϕ)) #Coriolis

    #Transport terms
    Γ_v = r*(ω^2)*cos(ϕ)*(sin(γ)*cos(ϕ)-cos(γ)*sin(ψ)*sin(ϕ))
    Γ_γ = (r/v)*(ω^2)*cos(ϕ)*(sin(γ)*sin(ψ)*sin(ϕ)+cos(γ)*cos(ϕ))
    Γ_ψ = -(r/v)*(ω^2)*(cos(ψ)*sin(ϕ)*cos(ϕ))/(cos(γ))

    #Attack angle profile. TO BE MODIFIED HERE.
    if  t <= 50.0
        α = 50.0
    else
        α = -t+110
    end

    #α = exp(-0.1*t)*cos(t) +15
    α = floor(Int, α)

    #Drag and Lift coeff/forces
    C_D = table_CD[α]
    C_L = table_CL[α]

    h = r - Re #m
    D = 0.5*exponential_atmosphere(h,ρ0,H)*v^2*(A_ref/m)*C_D #Drag Acceleration
    L = 0.5*exponential_atmosphere(h,ρ0,H)*v^2*(A_ref/m)*C_L

    #1-3 are same as simplified version
    ẋ[1] = v*sin(γ)
    ẋ[2] = (v/r)*(cos(γ)*cos(ψ)/cos(ϕ))
    ẋ[3] = (v/r)*cos(γ)*sin(ψ)
    ẋ[4] = -D+g_r*sin(γ)+g_ϕ*cos(γ)*sin(ψ) + Γ_v
    ẋ[5] = (1/v)*(L*cos(σ)+g_r*cos(γ)-g_ϕ*sin(γ)*sin(ψ))+(v/r)*cos(γ) + C_γ + Γ_γ
    ẋ[6] = -(1/(v*cos(γ)))*(L*sin(σ)-g_ϕ*cos(ψ)) -(v/r)*cos(γ)*cos(ψ)*tan(ϕ) +C_ψ + Γ_ψ
    return ẋ
end
