#full aero
function exponential_atmosphere(h)
    ρ0 = 0.0158 #kg/m3 here #*10^9 #0.026455*10^9 #sea level density (kg/km^3)
    h0 = 9.35458*1e3 #scale height on Mars(m)
    ρ = ρ0*exp(-h/h0)
end

a=1

#test for offline coefficients computation

function compute_aero2(δ, r_min, r_cone, r_G, v, ω, α, β)
    #general case: v is the relative velocity of spacecraft wrt atm
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    v_body = v*[sin(α)*cos(β), -sin(β), cos(α)*cos(β)] #minus on the y because of my body frame
    #@show(norm(v))
    #@show(norm(v_body))
    dr = 0.01
    dv = pi/20.0
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
            v_element = v_body - cross(ω, r_element-r_G)
            v_element = v_element/(norm(v_element))
            if n̂'*v_element > 0
                #dC = dot(n̂, v_body/norm(v_body))*dA*n̂
                #F_element = -0.5*exponential_atmosphere(h)*(norm(v_body)^2)*dC
                F_element = n̂'*v_element*dA*n̂*2*(n̂'*v_element) #CAREFUL
                τ_element = cross((r_element-r_G), F_element)
                CF_aero_body = CF_aero_body + F_element #m^2
                Cτ_aero_body = Cτ_aero_body + τ_element #m^3
            end
        end
    end
    return CF_aero_body/A_ref, Cτ_aero_body/(A_ref*L_ref)#CF_aero_body*(A_cone/A), Cτ_aero_body*(A_cone/A); #CF_aero_body, Cτ_aero_body  #
end

function table_aero(δ, r_min, r_cone, r_G, β)
    α = 0.0:1.0:359.0
    n = length(α)
    table_CF = zeros(n, 3)
    table_Cτ = zeros(n, 3)
    for i = 1:n
        CF, Cτ = compute_aero2(δ, r_min, r_cone, r_G, α[i]*pi/180, β)
        table_CF[i, :] = CF
        table_Cτ[i, :] = Cτ
    end
    return table_CF, table_Cτ
end

function compute_aero_sphere(δ, r_min, r_cone, r_G, v, ω, α, β)
    #general case: v is the relative velocity of spacecraft wrt atm
    CF_aero_body = [0.0;0.0;0.0]
    Cτ_aero_body = [0.0;0.0;0.0]
    v_body = v*[sin(α)*cos(β), -sin(β), cos(α)*cos(β)] #okay convention chosen, fair
    l = (r_cone-r_min)/tan(δ)
    dz = 0.01
    dv = pi/20.0
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
            v_element = v_body - cross(ω, r_element-r_G)
            v_element = v_element/(norm(v_element))
            if n̂'*v_element > 0
                F_element = n̂'*v_element*dA*n̂*2*(n̂'*v_element)
                τ_element = cross((r_element-r_G), F_element) #okay computed at COM
                CF_aero_body = CF_aero_body + F_element #m^2
                Cτ_aero_body = Cτ_aero_body + τ_element #m^3
            end
        end
    end
    return CF_aero_body/A_ref, Cτ_aero_body/(A_ref*L_ref)  #CF_aero_body*(A_sphere/A), Cτ_aero_body*(A_sphere/A)
end

function coeff_aero_spherecone(δ, r_min, r_cone, r_G, v, ω, α, β)
    #Based on traditional notation, here this function gives CA, -CY (=-CS), CN (or -CX, -CY, -CZ)
    #BUT as I changed the order, I actually return in this order: CN, -CY, CA (That's for FORCES)
    r_sphere = r_min/cos(δ) #radius where cone is cut
    A_ref_s = pi*r_sphere^2
    L_ref_s = r_sphere
    A_ref_c = pi*r_cone^2
    L_ref_c = r_cone
    A_ref = pi*r_cone^2
    L_ref = r_cone
    CF_cone, Cτ_cone = compute_aero2(δ, r_min, r_cone, r_G, v, ω, α, β)
    CF_sphere, Cτ_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, v, ω, α, β)
    CF = (CF_cone*A_ref_c+CF_sphere*A_ref_s)/A_ref
    Cτ = (Cτ_cone*A_ref_c*L_ref_c+Cτ_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
    return CF, Cτ
end

function table_aero_spherecone(δ, r_min, r_cone, r_G, β)
    #Based on traditional notation, here this function gives CA, -CY (=-CS), CN (or -CX, -CY, -CZ)
    #BUT as I changed the order, I actually return in this order: CN, -CY, CA (That's for FORCES)
    α = 0.0:1.0:359.0
    n = length(α)
    table_CF = zeros(n, 3)
    table_Cτ = zeros(n, 3)
    r_sphere = r_min/cos(δ) #radius where cone is cut
    A_ref_s = pi*r_sphere^2
    L_ref_s = r_sphere
    A_ref_c = pi*r_cone^2
    L_ref_c = r_cone
    A_ref = pi*r_cone^2
    L_ref = r_cone
    for i = 1:n
        CF_cone, Cτ_cone = compute_aero2(δ, r_min, r_cone, r_G, α[i]*pi/180, β)
        CF_sphere, Cτ_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, α[i]*pi/180, β)
        table_CF[i, :] = (CF_cone*A_ref_c+CF_sphere*A_ref_s)/A_ref
        table_Cτ[i, :] = (Cτ_cone*A_ref_c*L_ref_c+Cτ_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
    end
    return table_CF, table_Cτ
end

function drag_lift_coeff(δ, r_min, r_cone, r_G, α, β)
    r_sphere = r_min/cos(δ) #radius where cone is cut
    A_ref_s = pi*r_sphere^2
    L_ref_s = r_sphere
    A_ref_c = pi*r_cone^2
    L_ref_c = r_cone
    A_ref = pi*r_cone^2
    L_ref = r_cone
    CF_cone, Cτ_cone = compute_aero2(δ, r_min, r_cone, r_G, α, β)
    CF_sphere, Cτ_sphere = compute_aero_sphere(δ, r_min, r_cone, r_G, α, β)
    CF = (CF_cone*A_ref_c+CF_sphere*A_ref_s)/A_ref #this is CN, -CY, CA (in this order)
    Cτ = (Cτ_cone*A_ref_c*L_ref_c+Cτ_sphere*A_ref_s*L_ref_s)/(A_ref*L_ref)
    CN = CF[1]
    CA = CF[3]
    CD = CA*cos(α)+CN*sin(α)
    CL = CN*cos(α)-CA*sin(α)
    return CD, CL
end

function drag_lift_table(δ, r_min, r_cone, r_G, β)
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
β = 20.0*pi/180.0
δ = 70.0*pi/180
r_min = 0.09144*cos(δ)
r_cone = 0.762/2
r_G = [0.0;0.0;0.05]
t1, t2 = table_aero_spherecone(δ, r_min, r_cone, r_G, β)
α = 0.0:1.0:60.0
Plots.scatter(α, t2[1:61, 2])

#Forces coefficients
p1 = Plots.plot(α, t1[1:61, 1])
p2 = Plots.plot(α, t1[1:61, 2])
p3 = Plots.plot(α, t1[1:61, 3])
Plots.plot(p3, p2, p1, layout = (1, 3), legend = false)

#Moment coefficients
p1 = Plots.plot(α, -t2[1:61, 1])
p2 = Plots.plot(α, t2[1:61, 2])
p3 = Plots.plot(α, -t2[1:61, 3])
Plots.plot(p3, p2, p1, layout = (1, 3), legend = false)

##cone only
β = 20.0*pi/180.0
δ = 5.0*pi/180
r_min = 0.0
r_cone = 1.3
r_G = [0.0;0.0;0.0]
t1, t2 = table_aero(δ, r_min, r_cone, r_G, β)

#signs are such as to compare with reference in paper (y reversed and here normal outward, whereas paper normal inward
#and paper plots CA, CS, CN but Cl, Cm, Cn)
##force
α = 0.0:1.0:60.0
p1 = Plots.plot(α, t1[1:61, 1])
p2 = Plots.plot(α, t1[1:61, 2])
p3 = Plots.plot(α, t1[1:61, 3])
Plots.plot(p3, p2, p1, layout = (1, 3), legend = false)

##moment
α = 0.0:1.0:60.0
p1 = Plots.plot(α, -t2[1:61, 1])
p2 = Plots.plot(α, t2[1:61, 2])
p3 = Plots.plot(α, -t2[1:61, 3])
Plots.plot(p3, p2, p1, layout = (1, 3), legend = false) =#
