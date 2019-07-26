#Foot print plot on Mars
using MeshCat
using CoordinateTransformations
using GeometryTypes: GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh
using Colors: RGBA, RGB
using MeshIO
using FileIO
using LinearAlgebra
using ODE
using HCubature
using StaticArrays
using WebIO


Re = 3389.5 #3396.2 [km] Radius of the planet

function qmult(q1,q2)
    #scalar first
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

function compute_aero(δ, r_cone, r_G, r, v, q)
    #general case: v is the relative velocity of spacecraft wrt atm
    R_e = 3389.5
    h = norm(r)-R_e #altitude for atmosphere computation
    F_aero_body = [0.0;0.0;0.0]
    τ_aero_body = [0.0;0.0;0.0]
    v_body = qrot(qconj(q), v) #velocity of the spacecraft in body_frame in km.s-1
    #@show(v_body)
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
                τ_element = cross((r_element-r_G), F_element)*10^(3) #okay
                F_aero_body = F_aero_body + F_element
                τ_aero_body = τ_aero_body + τ_element
            end
        end
    end
    return F_aero_body*(A_cone/A), τ_aero_body*(A_cone/A)
end

function exponential_atmosphere(h)
    ρ0 = 0.0158*10^9 #0.026455*10^9 #sea level density (kg/km^3)
    h0 = 9.3545 #scale height on Mars(km)
    ρ = ρ0*exp(-h/h0)
end

#dynamics

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
            21.84 0.0 77.4]
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

    ω_mars = [0; 0; 7.095*10^(-5)] #[rad.s-1]
    v_rel = (v-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm

    #v_body = qrot(conj(q), v_rel)
    #α = atan(v_body[1], v_body[3])

    #control is a torque about the body z axis expressed in body frame at G global COM

    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q)
    F_aero_eci = qrot(q, F_aero_body)

    Cp = 1.0 #pitch damping coefficient
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

function dyn(t, x)
    return dyna(t, x, [0.0]) # torques input should be in km^2*kg*s^(-2) so should be small values
end

#=
function mat2quat(R)
    Q = zeros(4)
    r11, r12, r13 = R[1,1], R[1,2], R[1,3]
    r21, r22, r23 = R[2,1], R[2,2], R[2,3]
    r31, r32, r33 = R[3,1], R[3,2], R[3,3]
    N_0 = sqrt(1+r11+r22+r33)
    Q_0 = 0.5*[N_0; (r23-r32)/N_0; (r31-r13)/N_0; (r12-r21)/N_0]
    N_1 = sqrt(1+r11-r22-r33)
    Q_1 = 0.5*[(r23-r32)/N_1; N_1;(r12+r21)/N_1;(r31+r13)/N_1]
    N_2 = sqrt(1-r11+r22-r33)
    Q_2 = 0.5*[(r31-r13)/N_2;(r12+r21)/N_2; N_2;(r23+r32)/N_2]
    N_3 = sqrt(1-r11-r22+r33)
    Q_3 = 0.5*[(r12-r21)/N_3; (r31+r13)/N_3; (r23+r32)/N_3; N_3]
    if r22>=-r33 && r11>=-r22 && r11>=-r33
        Q = Q_0
    elseif r22<=-r33 && r11>=r22 && r11>=r33
        Q = Q_1
    elseif r22>=r33 && r11<=r22 && r11<=-r33
        Q = Q_2
    elseif r22<=r33 && r11<=-r22 && r11<=r33
        Q = Q_3
    end
    return Q
end =#


function mat2quat(M)
    #takes a rotation matrix as an input and returns a quaternion
    q1 = 0.5 * sqrt(1 + M[1,1] + M[2,2] + M[3,3])
    q2 = (M[2,3] - M[3, 2]) / (4 * q1)
    q3 = (M[3,1] - M[1,3]) / (4 * q1)
    q4 = (M[1,2] - M[2,1]) / (4 * q1)

    q = [q1, q2, q3, q4]
    q = q/norm(q)
    return q
end


#=
θ = 150*pi/180
M = [[0.0,0.0,-1.0], [-sin(θ), cos(θ), 0.0], [cos(θ),sin(θ),0.0]]
Q = mat2quat(M)
#x0 = [(3389.5+125)/Re, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.00001]
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]=#

function integration(fnc, x0)
    Re = 3389.5
    #x0_bb = [6453.448873708138/Re; -0564.603617093604/Re; 0.0; 1.0; 0.0;0.0;0.0; 0.683661419343188; 7.814285780472696; 0.0; 0.1;0.0;0.0] #orbit 5 degrees before
    #x0_b = [6477.113353992618/Re; -0113.058434141366/Re; 0.0; 1.0;0.0;0.0;0.0;0.136899033611736;7.842940382904379; 0.0;0.0;0.0;0.0] #orbit 1 degree before so should land before
    t_sim, z_sim = ode45(fnc, x0, 0:1:500, points=:specified)
    n = 13 #number of states
    Z = zeros(n,length(t_sim))

    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end


#Plot footprints on the surface of the Moon

function footprint()
    state_end = zeros(36, 13)
    pos_end = zeros(36, 3)
    list = 0:10:350
    for i in 1:1:length(list)
        θ = list[i]*pi/180
        if i != 10
        #=M = [0.0 0.0 -1.0;
            -sin(θ) cos(θ) 0.0;
             cos(θ) sin(θ) 0.0] #Z body axis towards the planet =#
        M = [-sin(θ) cos(θ) 0.0;
             0.0 0.0 1;
             cos(θ) sin(θ) 0.0]  #Z body axis in the direction of the initial velocity
        Q = mat2quat(M)
        Q = qconj(Q)
        @show(Q)
        x0 = [(3389.5+125)/Re, 0.0/Re, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 4.0, 0.0, 0.0, 0.0, 0.0]
        t_sim, Z= integration(dyn, x0)
        state_end[i, :] = Z[:, end]'
        pos_end[i, :] = Z[1:3, end]'
        @show(i)
    end
    end
    return pos_end, state_end
end

pos_end, state_end = footprint()
pos_end, state_end

##########################
###### Visualization #####
##########################
vis = Visualizer()
open(vis)

#Plot Mars in MeshCat
image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "Mars.png"))
texture = Texture(image=image)
material = MeshLambertMaterial(map=texture)
planet = HyperSphere(Point(0.,0,0), 0.1)
geometry = planet
setobject!(vis["planet"], geometry, material)
settransform!(vis["planet"], LinearMap(AngleAxis(pi/2, 1, 0, 0))) #rotate Planet

#Points (Footprint)
red_material = MeshPhongMaterial(color=RGBA(1, 0, 0, 1.0))
green_material = MeshPhongMaterial(color=RGBA(0, 1, 0, 1.0))


for i in 1:1:36
    if i != 10
        point = HyperSphere(Point(pos_end[i, 1]*10, pos_end[i, 2]*10, pos_end[i, 3]*0.1), 0.0002)
        setobject!(vis["point_$i"], point, red_material)
    end
end

#Departure point
pt_start = HyperSphere(Point(1.036878595663077*0.1, 0.0, 0.0), 0.0004)
setobject!(vis["pt_start"], pt_start, green_material)



#pos_end, state_end

#pos_end[30,:]
