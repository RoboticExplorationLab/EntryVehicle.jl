function dyna_coeffoff_COM_on_axis(t, x, u)
    ẋ = zeros(13)

    #@show(t)

    m_cone = 200.0 #kg
    m_p = 400.0 #kg
    m = m_p + m_cone #total mass in kg
    mu = 4.282837*1e13 #m^3.s-2  #42828.37 #km^3.s-2
    #r_G = [0.0;0.0;0.3] #COM in body frame of global spacecraft m
    δ = 70*pi/180 #radians
    r_cone = 1.325 #m
    h_cone = r_cone/tan(δ) #m
    Re = 3389.5*1e3 #m
    A_ref = pi*r_cone^2 #m2
    L_ref = r_cone #m

    #Inertia in kg*m^2
    #J = inertia(r_cone, h_cone) #J_total at global COM
    #=J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4] #check inertia matrix
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]=#

    #=No mass offset - kg*m^2
    J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193] =#

    J = [184.180 0.0 0.0;
        0.0 184.180 0.0;
        0.0 0.0 230.0]
    Jinv = inv(J)

    r = x[1:3]
    q = x[4:7]/norm(x[4:7])
    v = x[8:10]
    ω = x[11:13]

    if norm(r) > 3389.5*1e3

    #Compute aerodnyamic forces
    ω_mars = [0; 0; 7.095*10^(-5)] #rad.s-1
    #ω_mars = [0.0; 0.0; 0.0]
    v_rel = (v-cross(ω_mars, r)) #velocity of spacecraft wrt atm
    v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame
    α = acos(v_body[3]/norm(v_body)) #in fact total angle of attack, positif necessarily
    #=if α != 0.0 #&& abs(v_body[3]/(norm(v_body)*sin(α))) <= 1.0
            sin_β = -v_body[2]/(norm(v_body)*sin(α))
            if abs(v_body[1]/(norm(v_body)*sin(α))) > 1.0
                A = sign(v_body[1]/(norm(v_body)*sin(α)))*1.0
                if sin_β >=0.0
                    β = acos(A)
                else
                    β = 2*pi-acos(A)
                end
            else
                if sin_β >=0.0
                    β = acos(v_body[1]/(norm(v_body)*sin(α)))
                else
                    β = 2*pi-acos(v_body[1]/(norm(v_body)*sin(α)))
                end
            end
    else
            β = 0.0 #WRONG BUT NO WAY TO DETERMINE β in this case.
    end
    if t == 0.0
    @show(β)
end =#
    #α = floor(Int, α*180/pi)+1 #degrees
    #@show(v_body)
    #@show(α-1)

    ϕ_w = atan(v_body[2], v_body[1])

    q_v_b = [cos(ϕ_w/2); 0.0; 0.0; sin(ϕ_w/2)]

    α = α*180/pi #degrees

    #@show(α)
    #@show(β*180/pi)

    CF = [table_aero_chebyshev(α, 0.0, 181.0, C_FX);
        table_aero_chebyshev(α, 0.0, 181.0, C_FY);
        table_aero_chebyshev(α, 0.0, 181.0, C_FZ)]

    Cτ = [table_aero_chebyshev(α, 0.0, 181.0, C_τX);
        table_aero_chebyshev(α, 0.0, 181.0, C_τY);
        table_aero_chebyshev(α, 0.0, 181.0, C_τZ)] #given the plotted coeff curve but...

    Cd = [table_aero_chebyshev(α, 0.0, 181.0, DX);
        table_aero_chebyshev(α, 0.0, 181.0, DY);
        table_aero_chebyshev(α, 0.0, 181.0, DZ)]

    Cτ_final = [Cτ[1] - ω[1]*(L_ref/norm(v_body))*Cd[1]; Cτ[2] - ω[2]*(L_ref/norm(v_body))*Cd[2]; Cτ[3] - ω[3]*(L_ref/norm(v_body))*Cd[3]]

    h = (norm(r)-Re)
    F_aero_v = -0.5*exponential_atmosphere(h)*((norm(v_rel))^2)*CF*A_ref
    τ_aero_v = -0.5*exponential_atmosphere(h)*((norm(v_rel))^2)*Cτ_final*A_ref*L_ref
    F_aero_eci = qrot(qmult(q, q_v_b), F_aero_v)
    τ_aero_body = qrot(q_v_b, τ_aero_v)

    #controller stuff
    q_ref = Q
    kd = 100.0 #30.0
    kp = -100.0 #-20.0
    q_err = qmult(qconj(q), q_ref) #perfect measurements
    τ_c = -kd*(ω)-kp*q_err[2:4]

    #Compute gravitation + J2 acceleration
    J2 = 1960.45*(10^(-6))
    F_grav_eci = -m*mu*r/((norm(r)^3))

    F_J2_eci=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

    #Cp = 0.5 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci #+ F_J2_eci
    τ_total_body = τ_aero_body #+ τ_c #- Cp*[ω[1:2]; 0.0] #+ [0.0;0.0;u[1]] + τ_c#computed at COM #[0.0;0.0;0.0]
    #@show(τ_total_body)

    ẋ[1:3] = v
    ẋ[4:7] = 0.5*qmult(q, [0; ω])
    ẋ[8:10] = (F_total_eci/m)
    ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
else
    ẋ = zeros(13)
end
    return ẋ
end

##

function dyna_coeffon_COM_on_axis(t, x, u)
    ẋ = zeros(13)

    @show(t)

    m_cone = 152.0 #kg
    m_p = 400.0 #kg
    m = m_p + m_cone #total mass in kg
    mu = 4.282837*1e13 #m^3.s-2  #42828.37 #km^3.s-2
    r_G = [0.0;0.0;-0.189] #COM in body frame of global spacecraft m
    δ = 70*pi/180 #radians
    r_cone = 1.325 #m
    r_min = 0.2 #m
    h_cone = r_cone/tan(δ) #m
    Re = 3389.5*1e3 #m
    A_ref = pi*r_cone^2 #m2
    L_ref = r_cone #m

    #Inertia in kg*m^2
    #J = inertia(r_cone, h_cone) #J_total at global COM
    #=J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4] #check inertia matrix
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]=#

    #No mass offset - kg*m^2
    J = [184.180 0.0 0.0;
          0.0 184.180 0.0;
          0.0 0.0 230.0]

    #Jinv = [0.012952 0.0 0.0;
    #        0.0 0.012952 0.0;
    #        0.0 0.0  0.00986193]

    Jinv = inv(J)

    r = x[1:3]
    q = x[4:7]/norm(x[4:7])
    v = x[8:10]
    ω = x[11:13]

    if norm(r) > 3389.5*1e3

    #Compute aerodnyamic forces
    ω_mars = [0; 0; 7.095*10^(-5)] #rad.s-1
    #ω_mars = [0.0; 0.0; 0.0]
    v_rel = (v-cross(ω_mars, r)) #velocity of spacecraft wrt atm
    v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame

    β = asin(-v_body[2]/(norm(v_body)))
    α = atan(v_body[1], v_body[3])

    CF, Cτ = coeff_aero_spherecone(δ, r_min, r_cone, r_G, norm(v), ω, α, β)

    h = (norm(r)-Re)
    F_aero_v = -0.5*exponential_atmosphere(h)*((norm(v_rel))^2)*CF*A_ref
    τ_aero_v = -0.5*exponential_atmosphere(h)*((norm(v_rel))^2)*Cτ*A_ref*L_ref
    F_aero_eci = qrot(q, F_aero_v)
    τ_aero_body = τ_aero_v

    #controller stuff
    #q_ref = Q
    #kd = 100.0 #30.0
    #kp = -100.0 #-20.0
    #q_err = qmult(qconj(q), q_ref) #perfect measurements
    #τ_c = -kd*(ω)-kp*q_err[2:4]

    #Compute gravitation + J2 acceleration
    J2 = 1960.45*(10^(-6))
    F_grav_eci = -m*mu*r/((norm(r)^3))

    F_J2_eci=[J2*r[1]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[2]/norm(r)^7*(6*r[3]-1.5*(r[1]^2+r[2]^2));
     J2*r[3]/norm(r)^7*(3*r[3]-4.5*(r[1]^2+r[2]^2))]

    #Cp = 0.5 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci #+ F_J2_eci
    τ_total_body = τ_aero_body #-10.0*ω#+ τ_c #- Cp*[ω[1:2]; 0.0] #+ [0.0;0.0;u[1]] + τ_c#computed at COM #[0.0;0.0;0.0]
    #@show(τ_total_body)

    ẋ[1:3] = v
    ẋ[4:7] = 0.5*qmult(q, [0; ω])
    ẋ[8:10] = (F_total_eci/m)
    ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
else
    ẋ = zeros(13)
end
    return ẋ
end
