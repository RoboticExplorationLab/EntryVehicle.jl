module EntryVehicle

using LinearAlgebra

function vehicle_dynamics!(ẋ, x, u)

    #Spherical Earth
    Re = 6371.0 #Earth radius (km)

    #Cone with symmetry about x-axis
    #Need to add mass offset - velocity is 16° from cone centerline
    m_cone = 600.0 #Cone mass (kg)
    r_cone = 1.3 #Cone radius (m)
    θ_cone = 70*pi/180 #Cone half-angle
    h_cone = r_cone/tan(θ_cone) #Cone height (m)
    l_cone = 0.2 #Distance to cone COM from point (m)
    Cp = 0.1 #pitch damping coefficient

    #Inertia in kg*m^2
    Jxx = (3.0/10.0)*m_cone*r_cone^2;
    Jyy = (3.0/20.0)*m_cone*(r_cone^2 + 4*h_cone^2)
    J = [Jxx 0 0;
         0 Jyy 0;
         0 0 Jyy]
    Jinv = [1/Jxx 0 0;
            0 1/Jyy 0;
            0 0 1/Jyy]

    r = x[1:3] #position vector in ECI frame (km)
    q = x[4:7]/norm(x[4:7]) #body to ECI rotation
    v = x[8:10] #ECI velocity vector (rad/s)
    ω = x[11:13] #body frame angular velocity vector (km/s)

    #Calculate aerodynamic forces + torques in body frame
    h = norm(r) - Re #height above sea level (km)
    ρ = exponential_atmosphere(h) #atmospheric density (kg/m^3)
    vb = qrot(qconj(q),v);
    Fb = [0; 0; 0]
    tau = [0; 0; 0]
    dr = 0.1;
    dθ = pi/10;
    for ri in dr:dr:r_cone
        for θ = 0:dθ:(2*pi-dθ)
            nb = [cos(θ); sin(θ); -h_cone/r_cone];
            nb = nb/norm(nb);
            rc = ri*[cos(θ); sin(θ); h_cone/r_cone] - r_com;
            dA = ri*dθ*sqrt(dr^2 + (dr*h_cone/r_cone)^2);
            if nb'*vb > 0
                Fi = -0.5*ρ*dA*norm(vb)*nb*(nb'*vb);
                Fb = Fb + Fi;
                τ = τ + cross(rc, Fi);
            end
        end
    end

    #Add pitch damping
    τ = τ - Cp*[0; ω[2:3]]

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
    q3 = [q1[1]*q2[1] - q1[2:4]'*q2[2:4]';
          q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4], q2[2:4])]
end

function qrot(q,x)
    xrot = x + 2.0*cross(q[2:4], cross(q[2:4],x) + q[1]*x)
end

function qconj(q)
    qc = [q[1]; -q[2:4]]
end

end # module
