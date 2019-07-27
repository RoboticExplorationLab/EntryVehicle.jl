using LinearAlgebra
using Plots
using ODE
using HCubature
using StaticArrays
using PyPlot
pyplot()

#Dynamics on Mars
#NEED TO CHANGE ATM MODEL FOR MARS

Re = 3389.5

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
    dr = 0.1
    dv = pi/10.0
    A = 0.0
    A_cone = (pi*r_cone*sqrt(r_cone^2 + (r_cone/tan(δ))^2))*10^(-6) #km^2
    for r_i in 0:dr:r_cone
        for v in 0:dv:(2*pi-dv)
            n̂ = [cos(v)*cos(δ); sin(v)*cos(δ); sin(δ)] #normal outward
            dA = r_i*dv*dr/sin(δ) #dA in m^2
            dA = dA*10^(-6) #dA in km^2
            r_element = [r_i*cos(v); r_i*sin(v); (r_cone-r_i)/tan(δ)]
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

function exponential_atmosphere(h)
    ρ0 = 0.026455*10^9 #sea level density (kg/km^3)
    h0 = 10.8 #scale height on Mars(km)
    ρ = ρ0*exp(-h/h0)
end

#Plots.plot Results


function plot_traj(X)
    Plots.plot(X[1, :]*Re, X[2, :]*Re, label="Spacecraft Trajectory")
end

function plot_altitude(X, t_sim)
    Re = 3389.5
    R = zeros(length(t_sim))
    for i = 1:length(t_sim)
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-1)*Re
    end
    Plots.plot(t_sim, R, label="Spacecraft Altitude")
end

function plot_quaternions(X)
    Plots.plot(X[4, :])
    Plots.plot!(X[5, :])
    Plots.plot!(X[6, :])
    Plots.plot!(X[7, :])
end

function plot_ang_vel(X)
    Plots.plot(X[11, :])
    Plots.plot!(X[12, :])
    Plots.plot!(X[13, :])
end

function norm_quaternions(X, t_sim)
    N = zeros(length(t_sim))
    for i = 1:length(t_sim)
        N[i] = norm(X[4:7, i])
    end
    Plots.plot(t_sim, N)
end

function plot_vel(X, t_sim)
        V = zeros(length(t_sim))

        for i = 1:length(t_sim)
            V[i] = sqrt(X[8, i]^2+X[9, i]^2+X[10, i]^2)
        end
        Plots.plot(t_sim, V, label="Spacecraft Velocity")
end

function vector_z_body(X, t_sim)
    Z_b = zeros(length(t_sim), 3)
    for i=1:length(t_sim)
        q = X[4:7, i]'
        Z_b[i, :] = qrot(qconj(q), [0.0;0.0;1.0])'
    end
    return(Z_b)
end

#dynamics


function inertia(r_cone, h_cone)
    m_cone = 200.0
    m_point = 400.0
    m = m_cone + m_point  J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]
    Jzz = (3.0/10.0)*m_cone*r_cone^2;
    Jxx = (3.0/20.0)*m_cone*(r_cone^2 + 4*h_cone^2)
    J = [Jxx 0 0;
         0 Jxx 0;
         0 0 Jzz] #cone at its center of gravity
    G_cone = [0.0; 0.0; 0.118] #meters
    G_point = [0.0; 0.0; 0.118]
    #G_point = [0.3; 0.0; 0.391] #meters
    G = (m_cone*G_cone+m_point*G_point)/m
    x = G_point - G_cone
    J_t_Gp = J-m_cone*[x[2]^2+x[3]^2 -x[1]*x[2] -x[1]*x[3];
                      -x[1]*x[2] x[1]^2+x[3]^2 -x[2]*x[3];
                      -x[1]*x[3] -x[2]*x[3] x[1]^2+x[2]^2]
    y = G - G_point
    J_t_G = J_t_Gp - m*[y[2]^2+y[3]^2 -y[1]*y[2] y[1]*y[3];
                      -y[1]*y[2] y[1]^2+y[3]^2 -y[2]*y[3];
                      -y[1]*y[3] -y[2]*y[3] y[1]^2+y[2]^2]
    Jinv = inv(J_t_G)
    return J_t_G, Jinv
end


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
    J = [57.3336 0.0 10.92;
            0.0 33.3336 0.0;
            21.84 0.0 77.4]
    Jinv = [0.0184324 0.0 -0.00260054;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537]

    #No mass offset
    #=J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193]=#

    r = x[1:3]
    q = x[4:7]
    v = x[8:10]
    ω = x[11:13]

    if norm(r) > 1.0

    ω_mars = [0; 0; 0.0040612*pi/180]
    v_rel = (v-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm

    #v_body = qrot(conj(q), v_rel)
    #α = atan(v_body[1], v_body[3])

    #control is a torque about the body z axis expressed in body frame at G global COM

    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q)
    F_aero_eci = qrot(q, F_aero_body)

    Cp = 1.0 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci
    τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;u[1]] #computed at COM

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
end

#=function mat2quat(M)
    #takes a rotation matrix as an input and returns a quaternion
    q1 = 0.5 * sqrt(1 + M[1][1] + M[2][2] + M[3][3])
    q2 = (M[2][3] - M[3][2]) / (4 * q1)
    q3 = (M[3][1] - M[1][3]) / (4 * q1)
    q4 = (M[1][2] - M[2][1]) / (4 * q1)

    q = [q1, q2, q3, q4]
    q = q/norm(q)
    return q
end =#

θ = 100*pi/180
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
@show(Q)
#x0 = [(3389.5+125)/Re, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.00001]
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

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

t_sim, Z = integration(dyn, x0)

Z[1, 200]
Z[2, 200]
Z[3, 200]


plot_traj(Z)
xlabel!("X [km]")
ylabel!("Y [km]")
title!("Spacecraft Trajectory in MCI - XY projection")


plot_altitude(Z, t_sim)
xlabel!("t [s]")
ylabel!("altitude [km]")
title!("Spacecraft Altitude")


plot_vel(Z, t_sim)
xlabel!("t [s]")
ylabel!("Velocity [km/s]")
title!("Spacecraft Velocity")


plot_quaternions(Z)
xlabel!("t [s]")
ylabel!("Quaternions")
title!("Spacecraft Quaternions")


plot_ang_vel(Z)
norm_quaternions(Z, t_sim)
#Plots.plot(z_sim) #Okay now everything has a similar scale

#Plots.plots in 3D
#=using3D()
pygui(true)

fig = figure()
ax = fig[:gca](projection="3d")

n = 100
u = range(0.0,2*π,length=n);
v = range(0.0,π,length=n);
Re = 3389.5 #planet radius
x = Re*cos.(u) * sin.(v)';
y = Re*sin.(u) * sin.(v)';
z = Re*ones(n) * cos.(v)';

surf(x,y,z, zorder = 2)
ax.Plots.plot(Z[1, :]*Re, Z[2, :]*Re, Z[3, :]*Re, color= :yellow, linewidth=3, zorder=1)=#


#Plots.plot footprints on the surface of the Moon

function footprint()
    state_end = zeros(36, 13)
    pos_end = zeros(36, 3)
    list = 0:10:350
    for i in 1:1:length(list)
        θ = list[i]*pi/180
        if i != 19
        M = [[0.0,0.0,-1.0], [-sin(θ), cos(θ), 0.0], [cos(θ),sin(θ),0.0]]
        Q = mat2quat(M)
        x0 = [(3389.5+125)/Re, 0.0/Re, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.00001]
        t_sim, Z= integration(dyn, x0)
        state_end[i, :] = Z[:, end]'
        pos_end[i, :] = Z[1:3, end]'
        @show(i)
    end
    end
    return pos_end, state_end
end

#pos_end, state_end = footprint()
#scatter(pos_end[:, 1]*Re, pos_end[:, 2]*Re, pos_end[:, 3]*Re, color = "red")


#=
using PyPlots.plot
using3D()

pygui(true)

fig = figure()
ax = fig[:gca](projection="3d")

n = 100
u = range(0.0,2*π,length=n);
v = range(0.0,π,length=n);

Re = 3000
x = Re*cos.(u) * sin.(v)';
y = Re*sin.(u) * sin.(v)';
z = Re*ones(n) * cos.(v)';

surf(x, y, z)
ax.scatter(pos_end[:, 1]*Re, pos_end[:, 2]*Re, pos_end[:, 3]*Re, color = "red")=#
