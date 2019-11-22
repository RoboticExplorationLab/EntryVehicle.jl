#File contains functions commputing Aerodynamic Coefficients using Panel methods

####################################
####Pure Cone Geometry Procedure####
####################################

function compute_aero(δ, r_cone, r_G, r, v, q)
    #general case: v is the relative velocity of spacecraft wrt atm
    R_e = 3389.5
    h = norm(r)-R_e #altitude for atmosphere computation
    F_aero_body = [0.0;0.0;0.0]
    τ_aero_body = [0.0;0.0;0.0]
    v_body = qrot(qconj(q), v) #velocity of the spacecraft in body_frame in km.s-1
    #@show(norm(v))
    #@show(norm(v_body))
    #if norm(r) > 3489.0
    #@show(acos(v_body[3]/norm(v_body)))
    #end
    dr = 0.1
    dv = pi/10.0
    A = 0.0
    A_cone = (pi*r_cone*sqrt(r_cone^2 + (r_cone/tan(δ))^2))*10^(-6) #km^2
    for r_i in 0:dr:r_cone
        for a in 0:dv:(2*pi-dv)
            n̂ = [cos(a)*cos(δ); sin(a)*cos(δ); sin(δ)] #normal outward
            dA = r_i*dv*dr/sin(δ) #dA in m^2pos_end, state_end

            dA = dA*10^(-6) #dA in km^2
            r_element = [r_i*cos(a); r_i*sin(a); (r_cone-r_i)/tan(δ)]
            A = A + dA
            if n̂'*v_body > 0
                #dC = dot(n̂, v_body/norm(v_body))*dA*n̂
                #F_element = -0.5*exponential_atmosphere(h)*(norm(v_body)^2)*dC
                F_element = -0.5*exponential_atmosphere(h)*n̂'*v_body*dA*n̂*norm(v_body)*2*(n̂'*v_body/norm(v_body)) #CAREFUL
                τ_element = cross((r_element-r_G), F_element)*10^(3)
                F_aero_body = F_aero_body + F_element
                τ_aero_body = τ_aero_body + τ_element
            end
        end
    end
    return F_aero_body*(A_cone/A), τ_aero_body*(A_cone/A)
end

#This function is for a cut cone up to r_min
function compute_aero2(δ, r_min, r_cone, r_G, α)
    #general case: v is the relative velocity of spacecraft wrt atm
    R_e = 3389.5
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    v_body = [sin(α), 0.0, cos(α)]
    #@show(norm(v))
    #@show(norm(v_body))
    dr = 0.1
    dv = pi/10.0
    A = 0.0
    A_cone = (pi*r_cone*sqrt(r_cone^2 + (r_cone/tan(δ))^2))*10^(-6) #km^2
    for r_i in r_min:dr:r_cone
        for a in 0:dv:(2*pi-dv)
            n̂ = [cos(a)*cos(δ); sin(a)*cos(δ); sin(δ)] #normal outward
            dA = r_i*dv*dr/sin(δ) #dA in m^2pos_end, state_end
            dA = dA*10^(-6) #dA in km^2
            r_element = [r_i*cos(a); r_i*sin(a); (r_cone-r_i)/tan(δ)]
            A = A + dA
            if n̂'*v_body > 0
                #dC = dot(n̂, v_body/norm(v_body))*dA*n̂
                #F_element = -0.5*exponential_atmosphere(h)*(norm(v_body)^2)*dC
                F_element = n̂'*v_body*dA*n̂*2*(n̂'*v_body) #CAREFUL
                τ_element = cross((r_element-r_G), F_element)*10^(3)
                CF_aero_body = CF_aero_body + F_element #km
                Cτ_aero_body = Cτ_aero_body + τ_element #m
            end
        end
    end
    return CF_aero_body, Cτ_aero_body #CF_aero_body*(A_cone/A), Cτ_aero_body*(A_cone/A) #
end

function table_aero(δ, r_min, r_cone, r_G)
    α = 0.0:1.0:359.0
    n = length(α)
    table_CF = zeros(n, 3)
    table_Cτ = zeros(n, 3)
    for i = 1:n
        CF, Cτ = compute_aero2(δ, r_min, r_cone, r_G, α[i]*pi/180)
        table_CF[i, :] = CF
        table_Cτ[i, :] = Cτ
    end
    return table_CF, table_Cτ
end

####################################
####Pure Spherical Cap Procedure####
####################################

function compute_aero_sphere(δ, r_min, r_cone, r_G, α)
    #general case: v is the relative velocity of spacecraft wrt atm
    R_e = 3389.5
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    v_body = [sin(α), 0.0, cos(α)] #okay convention chosen, fair
    l = (r_cone-r_min)/tan(δ)
    dz = 0.1
    dv = pi/10.0
    A = 0.0
    r_sphere = r_min/cos(δ) #radius where cone is cut
    h = r_sphere*(1-sin(δ)) #height above the cut
    A_sphere = pi*(h^2+r_min^2)*(10^(-6)) #km^2 # overall area of the spherical cap
    for z in r_sphere-h:dz:r_sphere #check that
        for a in 0:dv:(2*pi-dv)
            n̂ = [1*cos(a)*sqrt(1-(z/r_sphere)^2); sin(a)*sqrt(1-(z/r_sphere)^2); z/r_sphere] #normal outward
            dA = (sqrt(r_sphere^2-z^2)*dv)*dz #dA in m^2
            dA = dA*10^(-6) #dA in km^2
            r_element = [sqrt(r_sphere^2-z^2)*cos(a); sqrt(r_sphere^2-z^2)*sin(a); l-r_sphere+h+z]
            A = A + dA
            if n̂'*v_body > 0
                F_element = n̂'*v_body*dA*n̂*2*(n̂'*v_body) #CAREFUL
                τ_element = cross((r_element-r_G), F_element)*10^(3)
                CF_aero_body = CF_aero_body + F_element #km
                Cτ_aero_body = Cτ_aero_body + τ_element #m
            end
        end
    end
    return CF_aero_body, Cτ_aero_body  #CF_aero_body*(A_sphere/A), Cτ_aero_body*(A_sphere/A)
end

####################################
########Sphere Cone tables##########
####################################

function table_aero_spherecone(δ, r_min, r_cone, r_G)
    α = 0.0:1.0:359.0
    n = length(α)
    table_CF = zeros(n, 3)
    table_Cτ = zeros(n, 3)
    for i = 1:n
        CF_cone, Cτ_cone = compute_aero2(δ, r_min, r_cone, r_G, α[i]*pi/180)
        CF_sphere, Cτ_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, α[i]*pi/180)
        table_CF[i, :] = CF_cone+CF_sphere
        table_Cτ[i, :] = Cτ_cone+Cτ_sphere
    end
    return table_CF, table_Cτ
end

####################################
########Exp Atmospheric Model#######
####################################

function exponential_atmosphere(h)
    ρ0 = 0.0158*10^9 #0.026455*10^9 #sea level density (kg/km^3)
    h0 = 9.3545 #scale height on Mars(km)
    ρ = ρ0*exp(-h/h0)
end

#=test space
δ = 70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]

t1, t2 = table_aero(δ, r_min, r_cone, r_G)
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

α = 0.0:1.0:359.0
plot(α, table_CF[:, 2]) =#
