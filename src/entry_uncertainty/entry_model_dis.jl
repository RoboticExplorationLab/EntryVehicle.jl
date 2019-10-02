#File contains entry vehicle model for Mars EDL

function dyna_dis_temp(t, dt, x1, u, w)
    x2 = zeros(15)

    #@show(t)

    m_cone = 200.0 #kg
    m_p = 400.0 #kg
    m = m_p + m_cone #total mass in kg
    mu = 42828.37 #km^3.s-2
    r_G = [0.2;0.0;0.3] #COM in body frame of global spacecraft
    δ = 70*pi/180 #radians
    r_cone = 1.3 #in meters
    h_cone = r_cone/tan(δ)
    Re = 3389.5 #Radius in km

    #Inertia in kg*m^2
    #J = inertia(r_cone, h_cone) #J_total at global COM
    J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4] #check inertia matrix
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]

    #No mass offset
    #=J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193] =#

    r = x1[1:3]
    q = x1[4:7]/norm(x1[4:7])
    v = x1[8:10]
    ω = x1[11:13]

    if norm(r) > 1.0

    #Compute aerodnyamic forces
    ω_mars = [0; 0; 7.095*10^(-5)]
    v_rel = (v-cross(ω_mars, r*Re)) + w[2:4] #velocity of spacecraft wrt atm adding wind w

    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q, w[1])
    F_aero_eci = qrot(q, F_aero_body)

    #Compute gravitation + J2 acceleration
    x = r[1]
    y = r[2]
    z = r[3]
    J2 = 1960.45*(10^(-6))
    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    #=F_J2_X = (5*((z/norm(r))^2)-1)*x
    F_J2_Y = (5*((z/norm(r))^2)-1)*y
    F_J2_Z = (5*((z/norm(r))^2)-3)*z
    F_J2_eci = ((3*J2*mu*Re^2)/(2*(norm(r)^5)))*[F_J2_X; F_J2_Y; F_J2_Z] =#
    F_J2_eci=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

    Cp = 0.5 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci #+ F_J2_eci
    τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;u[1]] #computed at COM #[0.0;0.0;0.0]

    x2[1:3] = x1[1:3] + (v/Re)*dt
    x2[4:7] = x1[4:7] + 0.5*qmult(q, [0; ω])*dt
    x2[4:7] = x2[4:7]/(norm(x2[4:7]))
    x2[8:10] = x1[8:10] + (F_total_eci/m)*dt
    x2[11:13] = x1[11:13] + Jinv*(τ_total_body-cross(ω, J*ω))*dt

    #Temperature computation
    r_nose = 0.5
    ϵ = 0.8
    σ = 5.67*10^(-8) #stefan-Boltzmann constant
    k = 1.9027*10^(-4)
    ρ = exponential_atmosphere(w[1], (norm(x2[1:3])-1.0)*Re)
    q_s = k*sqrt(ρ*10^(-9)/r_nose)*norm(x2[8:10]*10^3)^3.04
    x2[14] = (q_s/(ϵ*σ))^0.25
    x2[15] = (norm(x2[8:10]-x1[8:10])/dt)*(10^(3))/(9.81) #normalized acceleration of the vehicle
else
    x2 = x1
end
    return x2
end


function dyna_dis(t, dt, x1, u, w)
    x2 = zeros(13)

    #@show(t)

    m_cone = 200.0 #kg
    m_p = 400.0 #kg
    m = m_p + m_cone #total mass in kg
    mu = 42828.37 #km^3.s-2
    r_G = [0.2;0.0;0.3] #COM in body frame of global spacecraft
    δ = 70*pi/180 #radians
    r_cone = 1.3 #in meters
    h_cone = r_cone/tan(δ)
    Re = 3389.5 #Radius in km

    #Inertia in kg*m^2
    #J = inertia(r_cone, h_cone) #J_total at global COM
    J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4] #check inertia matrix
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]

    #No mass offset
    #=J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193] =#

    r = x1[1:3]
    q = x1[4:7]/norm(x1[4:7])
    v = x1[8:10]
    ω = x1[11:13]

    if norm(r) > 1.0

    #Compute aerodnyamic forces
    ω_mars = [0; 0; 7.095*10^(-5)]
    v_rel = (v-cross(ω_mars, r*Re)) + w[2:4] #velocity of spacecraft wrt atm adding wind w

    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q, w[1])
    F_aero_eci = qrot(q, F_aero_body)

    #Compute gravitation + J2 acceleration
    x = r[1]
    y = r[2]
    z = r[3]
    J2 = 1960.45*(10^(-6))
    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    #=F_J2_X = (5*((z/norm(r))^2)-1)*x
    F_J2_Y = (5*((z/norm(r))^2)-1)*y
    F_J2_Z = (5*((z/norm(r))^2)-3)*z
    F_J2_eci = ((3*J2*mu*Re^2)/(2*(norm(r)^5)))*[F_J2_X; F_J2_Y; F_J2_Z] =#
    F_J2_eci=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

    Cp = 0.5 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci #+ F_J2_eci
    τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;u[1]] #computed at COM #[0.0;0.0;0.0]

    x2[1:3] = x1[1:3] + (v/Re)*dt
    x2[4:7] = x1[4:7] + 0.5*qmult(q, [0; ω])*dt
    x2[8:10] = x1[8:10] + (F_total_eci/m)*dt
    x2[11:13] = x1[11:13] + Jinv*(τ_total_body-cross(ω, J*ω))*dt

    #Temperature computation
    #r_nose = 0.5
    #ϵ = 0.8
    #σ = 5.67*10^(-8) #stefan-Boltzmann constant
    #k = 1.9027*10^(-4)
    #ρ = exponential_atmosphere(w[1], (norm(x2[1:3])-1.0)*Re)
    #q_s = k*sqrt(ρ*10^(-9)/r_nose)*norm(x2[8:10]*10^3)^3.04
    #x2[14] = (q_s/(ϵ*σ))^0.25
else
    x2 = x1
end
    return x2
end

function dyna_coeffoff(t, x)
    ẋ = zeros(13)

    #@show(t)

    m_cone = 200.0 #kg
    m_p = 400.0 #kg
    m = m_p + m_cone #total mass in kg
    mu = 42828.37 #km^3.s-2
    r_G = [0.2;0.0;0.3] #COM in body frame of global spacecraft
    δ = 70*pi/180 #radians
    r_cone = 1.3 #in meters
    h_cone = r_cone/tan(δ)
    Re = 3389.5 #Radius in km

    #Inertia in kg*m^2
    #J = inertia(r_cone, h_cone) #J_total at global COM
    J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4] #check inertia matrix
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]

    #No mass offset
    #=J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193] =#

    r = x[1:3]
    q = x[4:7]/norm(x[4:7])
    v = x[8:10]
    ω = x[11:13]

    if norm(r) > 1.0

    #Compute aerodnyamic forces
    ω_mars = [0; 0; 7.095*10^(-5)]
    v_rel = (v-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm
    v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame
    α = acos(v_body[3]/norm(v_body)) #radians
    α = floor(Int, α*180/pi)+1 #degrees

    CF = table_CF[α, :]
    Cτ = table_Cτ[α, :]
    h = (norm(r)-1.0)*Re
    F_aero_body = -0.5*exponential_atmosphere(h)*norm(v_rel)^2*CF
    τ_aero_body = -0.5*exponential_atmosphere(h)*norm(v_rel)^2*Cτ
    F_aero_eci = qrot(q, F_aero_body)

    #Compute gravitation + J2 acceleration
    x = r[1]
    y = r[2]
    z = r[3]
    J2 = 1960.45*(10^(-6))
    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    #=F_J2_X = (5*((z/norm(r))^2)-1)*x
    F_J2_Y = (5*((z/norm(r))^2)-1)*y
    F_J2_Z = (5*((z/norm(r))^2)-3)*z
    F_J2_eci = ((3*J2*mu*Re^2)/(2*(norm(r)^5)))*[F_J2_X; F_J2_Y; F_J2_Z] =#
    F_J2_eci=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

    Cp = 0.1 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci + F_J2_eci
    τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;0.0] #computed at COM #[0.0;0.0;0.0]

    ẋ[1:3] = v/Re
    ẋ[4:7] = 0.5*qmult(q, [0; ω])
    ẋ[8:10] = F_total_eci/m
    ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
else
    ẋ = zeros(13)
end
    return ẋ
end
