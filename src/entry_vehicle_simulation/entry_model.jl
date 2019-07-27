#File contains entry vehicle model for Mars EDL

function dyna(t, x, u)
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

    ω_mars = [0; 0; 7.095*10^(-5)]
    v_rel = (v-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm

    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q)
    F_aero_eci = qrot(q, F_aero_body)

    Cp = 0.5 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci
    τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;u[1]] #computed at COM #[0.0;0.0;0.0]

    ẋ[1:3] = v/Re
    ẋ[4:7] = 0.5*qmult(q, [0; ω])
    ẋ[8:10] = F_total_eci/m
    ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
else
    ẋ = zeros(13)
end
    return ẋ
end
