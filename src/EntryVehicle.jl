module EntryVehicle

using LinearAlgebra

function vehicle_dynamics!(ẋ, x, u, t)

    #Flat Earth
    Re = 6371.0 #Earth radius (km)
    g = 9.81 #(m/s^2)

    #Cone with symmetry about z-axis
    #Need to add mass offset - velocity is 16° from cone centerline
    m_cone = 600.0 #Cone mass (kg)
    r_cone = 1.3 #Cone radius (m)
    θ_cone = 70*pi/180 #Cone half-angle
    h_cone = r_cone/tan(θ_cone) #Cone height (m)
    A_cone = pi*r_cone*sqrt(r_cone^2 + h_cone^2) #(m^2)
    com_cone = [0; 0; 0.2] #Vector to cone COM from point (m)
    Cp = 1.0 #pitch damping coefficient

    #Inertia in kg*m^2
    Jzz = (3.0/10.0)*m_cone*r_cone^2;
    Jxx = (3.0/20.0)*m_cone*(r_cone^2 + 4*h_cone^2)
    J = [Jxx 0 0;
         0 Jxx 0;
         0 0 Jzz]
    Jinv = [1/Jxx 0 0;
            0 1/Jxx 0;
            0 0 1/Jzz]

    r = x[1:3] #position vector in ECI frame (km)
    q = x[4:7]/norm(x[4:7]) #body to ECI rotation
    v = x[8:10] #ECI velocity vector (km/s)
    ω = x[11:13] #body frame angular velocity vector (rad/s)

    #Calculate aerodynamic forces + torques in body frame
    h = norm(r) - Re #height above sea level (km)
    ρ = exponential_atmosphere(h) #atmospheric density (kg/m^3)
    vb = qrot(qconj(q),v);
    Fb = [0; 0; 0]
    τ = [0; 0; 0]
    A = 0;
    dr = 0.1;
    dθ = pi/10;
    for ri in dr:dr:r_cone
        for θ = 0:dθ:(2*pi-dθ)
            nb = [cos(θ); sin(θ); -h_cone/r_cone];
            nb = nb/norm(nb);
            rc = ri*[cos(θ); sin(θ); h_cone/r_cone] - com_cone;
            dA = ri*dθ*sqrt(dr^2 + (dr*h_cone/r_cone)^2);
            A = A + dA
            if nb'*vb > 0
                Fi = -0.5*ρ*dA*norm(vb)*nb*(nb'*vb);
                Fb = Fb + Fi;
                τ = τ + cross(rc, Fi);
            end
        end
    end
    Feci = qrot(q,(A_cone/A)*Fb) - [0; 0; g]

    #Add pitch damping
    τ = (A_cone/A)*τ - Cp*[ω[1:2]; 0]

    ẋ[1:3] = v
    ẋ[4:7] = 0.5*qmult(q,[0; ω])
    ẋ[8:10] = Feci/m_cone
    ẋ[11:13] = Jinv*(τ - cross(ω, J*ω))

end

function exponential_atmosphere(h)
    ρ0 = 1.2 #sea level density (kg/m^3)
    h0 = 8.0 #scale height (km)
    ρ = ρ0*exp(-h/h0)
end

function spherical_gravity(r)
    mu = 398600.0 #standard gravitational parameter (km^3/s^2)
    g = -mu*r/(norm(r)^3)
end

function qmult(q1,q2)
    s1 = q1[1]
    s2 = q2[1]
    v1 = q1[2:4]
    v2 = q2[2:4]
    q3 = [s1*s2 - v1'*v2;
          s1*v2 + s2*v1 + cross(v1, v2)]
end

function qrot(q,x)
    xrot = x + 2.0*cross(q[2:4], cross(q[2:4],x) + q[1]*x)
end

function qconj(q)
    qc = [q[1]; -q[2:4]]
end

end # module
