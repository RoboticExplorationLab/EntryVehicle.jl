#Space Mechanics functions for orbit computation and stuff like that
##WARNINGG : need mu in meters and Re for some of the dynamical models

function ecef2geoc(r_ecef)
    #goes from the planet centered reference frame to the geocentric coordinates
    #takes the r_ecef vector with no scaling, in km
    Re = 3389.5;
    r_x = r_ecef[1]
    r_y = r_ecef[2]
    r_z = r_ecef[3]
    r = norm(r_ecef)
    h = r-Re
    ϕ = rad2deg(asin(r_z/r)) #latitude (from equator)
    λ = rad2deg(atan(r_y, r_x)) #longitude (from ECEF x axis)
    return [h; ϕ; λ]
end

a=1

function geoc2geod(r_geoc)
    #r_geoc = [h; ϕ; λ]
    Re = 3389.5
    e_g = 5.889*10^(-3)
    r_ecef = geoc2ecef(r_geoc)
    λ = r_geoc[3]
    ϕ = r_geoc[2]
    r_x = r_ecef[1]
    r_y = r_ecef[2]
    r_z = r_ecef[3]
    λ_p = λ
    ϕ_p = ϕ
    δ = 1
    ϵ = 10^(-10)
    r_xy = sqrt(r_x^2+r_y^2)
    while abs(δ) > ϵ
        N = Re/sqrt(1-(e_g^2)*sind(ϕ_p)^2)
        ϕ_p_new = atand((r_z+N*(e_g^2)*sind(ϕ_p)/r_xy))
        δ = ϕ_p_new-ϕ_p
        ϕ_p = ϕ_p_new

    end
    h_p = (r_xy/cosd(ϕ_p))-N
    return [h_p; ϕ_p; λ_p] #r_geod
end

function ecef2geod(r_ecef)
    r_geoc = ecef2geoc(r_ecef)
    r_geod = geoc2geod(r_geoc) #degree
    return r_geod
end


function eci2ecef(θ, r_eci)
    # θ is the GMST, given in radians here
    R = [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0] #as theta evolves with time, the matrix will be time dependent
    r_ecef = R*r_eci
    return r_ecef
end

a = 1

function oe2eci(r_oe)
    mu = 42828.37 #km^3.s-2
    #input angles in degrees
    a = r_oe[1]
    e = r_oe[2]
    i = r_oe[3]
    W = r_oe[4]
    w = r_oe[5]
    nu = r_oe[6]
    R11 = cosd(W)*cosd(w)-sind(W)*sind(w)*cosd(i);
    R12 = -cosd(W)*sind(w)-sind(W)*cosd(w)*cosd(i);
    R13 = sind(W)*sind(i);
    R21 = sind(W)*cosd(w)+cosd(W)*sind(w)*cosd(i);
    R22 = -sind(W)*sind(w)+cosd(W)*cosd(w)*cosd(i);
    R23 = -cosd(W)*sind(i);
    R31 = sind(w)*sind(i);
    R32 = cosd(w)*sind(i);
    R33 = cosd(i);
    R = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
    p = a*(1-e^2);
    r = p/(1+e*cosd(nu));
    x = r*cosd(nu);
    y = r*sind(nu);
    r_PQW = [x;y;0.0];
    v_PQW = [-sqrt(mu/p)*sind(nu);sqrt(mu/p)*(e+cosd(nu));0.0]
    r_ECI = R*r_PQW;
    v_ECI = R*v_PQW;
    return [r_ECI; v_ECI]
end

a=1

function eci2oe(r_eci, v_eci)
    mu = 42828.37 #398600.4418 42828.37; #[km3.s-2]
    v_norm = norm(v_eci);
    r_norm = norm(r_eci);
    h = cross(r_eci, v_eci); #vector h here
    h_norm = norm(h);
    energy = ((v_norm^2)/2)-(mu/r_norm);
    a = -mu/(2*energy);
    p = (h_norm^2)/mu;
    e_vec = (1/mu)*(((v_norm^2)-mu/r_norm)*r_eci-dot(r_eci,v_eci)*v_eci); #e_vec
    e = norm(e_vec);
    n = cross([0.0;0.0;1.0], h);
    i = acosd(h[3]/h_norm); #i is always less than 180°
    if e==0.0 && i==0.0 || e==0.0 && i==180.0
        w = NaN;
        W = NaN;
        l = acosd(r_eci[1]/r_norm);
        nu_0 = l;
    elseif e==0.0 && i!=0.0 || e==0.0 && i!=180.0
        w = NaN;
        if n[2]>0.0
        W = acosd(n[1]/norm(n));
        else
        W = 360.0-acosd(n[1]/norm(n));
        end
        if r[3]>0.0
            u = acosd((dot(n,r_eci))/(norm(n)*r_norm));
        else
            u = 360.0-acosd((dot(n,r_eci))/(norm(n)*r_norm));
        nu_0 = u;
        end
    elseif e!=0.0 && i==0.0 || e!=0.0 && i==180.0
        W = NaN;
        w = NaN;
        if dot(r_eci,v_eci)>0.0
            nu_0 = acosd((dot(e_vec, r_eci))/(e*r_norm));
        else
            nu_0 = 360.0 - acosd((dot(e_vec, r_eci))/(e*r_norm));
        end
    else
        if n[2]>0.0
            W = acosd(n[1]/norm(n));
        else
            W = 360.0-acosd(n[1]/norm(n));
        end
        if e_vec[3]>0.0
            w = acosd((dot(n,e_vec))/(norm(n)*e));
        else
            w = 360.0 - acosd((dot(n,e_vec))/(norm(n)*e));
        end
        if dot(r_eci,v_eci)>0.0
            nu_0 = acosd((dot(e_vec, r_eci))/(e*r_norm));
        else
            nu_0 = 360.0 - acosd((dot(e_vec, r_eci))/(e*r_norm));
        end
    end
    return [a;e;i;W;w;nu_0]
end

function M2E(M, e) #E, e of the orbit
    E = M #init
    δ = 1
    while abs(δ)>10^(-10)
        temp = E
        E = E-((E-e*sin(E)-M)/(1-e*cos(E))) #radians
        δ = E-temp
    end
    return E
end

#test space
#=using LinearAlgebra
r_ECI_0 = [-2827.4;-8916.9; -2340.8]; #[km]
v_ECI_0 = [4.6567; 3.2251; -5.9629]; #[km.s-1]
X = eci2oe(r_ECI_0, v_ECI_0)=#
