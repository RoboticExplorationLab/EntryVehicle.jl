
function exponential_atmosphere(h)
    ρ0 = 0.0158 #kg/m3 here #*10^9 #0.026455*10^9 #sea level density (kg/km^3)
    h0 = 9.35458*1e3 #scale height on Mars(m)
    ρ = ρ0*exp(-h/h0)
end

function speed_sound(h)
    #returns speed of sound with respect to altitude (fit model for h in m)
    #h in M
    #martian atmopshere data
    β0 = 223.8
    β1 = -0.0002004
    β2 = -1.588e-8
    β3 = 1.404e-13
    v_s = β0 + β1*h + β2*h^2 + β3*h^3
    return v_s
end

a=1

#test for offline coefficients computation

function compute_aero2(δ, r_min, r_cone, r_G, α)
    #general case: v is the relative velocity of spacecraft wrt atm
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    Cτ_damping = [0.0;0.0;0.0]
    v_body = [sin(α), 0.0, cos(α)]
    #@show(norm(v))
    #@show(norm(v_body))
    dr = 0.1
    dv = pi/10.0
    A = 0.0
    #A_cone = (pi*r_cone*sqrt(r_cone^2 + (r_cone/tan(δ))^2)) #m^2
    A_ref = pi*r_cone^2
    L_ref = r_cone
    for r_i in r_min:dr:r_cone
        for a in 0:dv:(2*pi-dv)
            n̂ = [cos(a)*cos(δ); sin(a)*cos(δ); sin(δ)] #normal outward
            dA = r_i*dv*dr/sin(δ) #dA in m^2pos_end, state_end
            dA = dA #dA in m^2
            r_element = [r_i*cos(a); r_i*sin(a); (r_cone-r_i)/tan(δ)]
            A = A + dA
            if n̂'*v_body > 0
                #dC = dot(n̂, v_body/norm(v_body))*dA*n̂
                #F_element = -0.5*exponential_atmosphere(h)*(norm(v_body)^2)*dC
                F_element = n̂'*v_body*dA*n̂*2*(n̂'*v_body) #CAREFUL
                τ_element = cross((r_element-r_G), F_element)
                CF_aero_body = CF_aero_body + F_element #m^2
                Cτ_aero_body = Cτ_aero_body + τ_element #m^3
                Cτ_damping[1] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[1])^2)*dA)/(A_ref*(L_ref^2))
                Cτ_damping[2] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[2])^2)*dA)/(A_ref*(L_ref^2))
                Cτ_damping[3] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[3])^2)*dA)/(A_ref*(L_ref^2))
            end
        end
    end
    return CF_aero_body/A_ref, Cτ_aero_body/(A_ref*L_ref), Cτ_damping#CF_aero_body*(A_cone/A), Cτ_aero_body*(A_cone/A); #CF_aero_body, Cτ_aero_body  #
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

function compute_aero_sphere(δ, r_min, r_cone, r_G, α)
    #general case: v is the relative velocity of spacecraft wrt atm
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    Cτ_damping = [0.0;0.0;0.0]
    v_body = [sin(α), 0.0, cos(α)] #okay convention chosen, fair
    l = (r_cone-r_min)/tan(δ)
    dz = 0.1
    dv = pi/10.0
    A = 0.0
    r_sphere = r_min/cos(δ) #radius where cone is cut
    h = r_sphere*(1-sin(δ)) #height above the cut
    A_sphere = pi*(h^2+r_min^2) #m^2 # overall area of the spherical cap
    A_ref = pi*r_sphere^2
    L_ref = r_sphere
    for z in r_sphere-h:dz:r_sphere #check that
        for a in 0:dv:(2*pi-dv)
            n̂ = [1*cos(a)*sqrt(1-(z/r_sphere)^2); sin(a)*sqrt(1-(z/r_sphere)^2); z/r_sphere] #normal outward
            dA = (sqrt(r_sphere^2-z^2)*dv)*dz #dA in m^2
            dA = dA #dA in m^2
            r_element = [sqrt(r_sphere^2-z^2)*cos(a); sqrt(r_sphere^2-z^2)*sin(a); l-r_sphere+h+z]
            A = A + dA
            if n̂'*v_body > 0
                F_element = n̂'*v_body*dA*n̂*2*(n̂'*v_body)
                τ_element = cross((r_element-r_G), F_element) #okay computed at COM
                CF_aero_body = CF_aero_body + F_element #m^2
                Cτ_aero_body = Cτ_aero_body + τ_element #m^3
                Cτ_damping[1] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[1])^2)*dA)/(A_ref*(L_ref^2))
                Cτ_damping[2] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[2])^2)*dA)/(A_ref*(L_ref^2))
                Cτ_damping[3] += (-4.0*(n̂'*v_body)*((cross(r_element-r_G, n̂)[3])^2)*dA)/(A_ref*(L_ref^2))
            end
        end
    end
    return CF_aero_body/A_ref, Cτ_aero_body/(A_ref*L_ref), Cτ_damping #CF_aero_body*(A_sphere/A), Cτ_aero_body*(A_sphere/A)
end

function table_aero_spherecone(δ, r_min, r_cone, r_G)
    #Based on traditional notation, here this function gives CA, -CY, CN (or -CX, -CY, -CZ)
    #BUT as I changed the order, I actually return in this order: CN, -CY, CA (That's for FORCES)
    α = 0.0:1.0:359.0
    n = length(α)
    table_CF = zeros(n, 3)
    table_Cτ = zeros(n, 3)
    table_damping = zeros(n, 3)
    r_sphere = r_min/cos(δ) #radius where cone is cut
    A_ref_s = pi*r_sphere^2
    L_ref_s = r_sphere
    A_ref_c = pi*r_cone^2
    L_ref_c = r_cone
    A_ref = pi*r_cone^2
    L_ref = r_cone
    for i = 1:n
        CF_cone, Cτ_cone, Cd_cone = compute_aero2(δ, r_min, r_cone, r_G, α[i]*pi/180)
        CF_sphere, Cτ_sphere, Cd_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, α[i]*pi/180)
        table_CF[i, :] = (CF_cone*A_ref_c+CF_sphere*A_ref_s)/A_ref
        table_Cτ[i, :] = (Cτ_cone*A_ref_c*L_ref_c+Cτ_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
        table_damping[i, :] = (Cd_cone*A_ref_c*L_ref_c+Cd_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
    end
    return table_CF, table_Cτ, table_damping
end

function drag_lift_coeff(δ, r_min, r_cone, r_G, α)
    r_sphere = r_min/cos(δ) #radius where cone is cut
    A_ref_s = pi*r_sphere^2
    L_ref_s = r_sphere
    A_ref_c = pi*r_cone^2
    L_ref_c = r_cone
    A_ref = pi*r_cone^2
    L_ref = r_cone
    CF_cone, Cτ_cone = compute_aero2(δ, r_min, r_cone, r_G, α)
    CF_sphere, Cτ_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, α)
    CF = (CF_cone*A_ref_c+CF_sphere*A_ref_s)/A_ref #this is CN, -CY, CA (in this order)
    Cτ = (Cτ_cone*A_ref_c*L_ref_c+Cτ_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
    CN = CF[1]
    CA = CF[3]
    CD = CA*cos(α)+CN*sin(α)
    CL = CN*cos(α)-CA*sin(α)
    return CD, CL
end

function drag_lift_table(δ, r_min, r_cone, r_G)
    α = 0.0:1.0:180.0
    n = length(α)
    table_CD = zeros(n)
    table_CL = zeros(n)
    for i=1:length(α)
        CD, CL = drag_lift_coeff(δ, r_min, r_cone, r_G, α[i]*pi/180)
        table_CD[i] = CD
        table_CL[i] = CL
    end
    return table_CD, table_CL
end

#=test space
δ = 70*pi/180
r_min = 0.09144*cos(δ)
r_cone = 0.762/2
r_G = [0.0; 0.0; -0.1]

t1, t2 = table_aero(δ, r_min, r_cone, r_G)
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

td, tl = drag_lift_table(δ, r_min, r_cone, r_G)
Plots.plot(tl)

#Sphere cone in the case of an axisymmetric body makes sense (same curve as the hypersonic paper)
#Difference certainly due to the references chosen and their Beta angle that I did not consider (because symmetric totally)
t1, t2 = table_aero_spherecone(δ, r_min, r_cone, r_G)
α = 0.0:1.0:60.0
Plots.scatter(α, t2[1:61, 2])

p1 = Plots.plot(α, t1[1:61, 1])
Plots.xlabel!("α")
p2 = Plots.plot(α, t1[1:61, 3])
p3 = Plots.plot(α, t2[1:61, 2])
Plots.plot(p1, p2, p3, layout = (1, 3), legend = false)
savefig("aero_coeff") =#

#=speed of sound tests
H = 0.0:1000.0:150000.0
V = [speed_sound(h) for h in H]
Plots.plot(H/(1e3), V) =#
