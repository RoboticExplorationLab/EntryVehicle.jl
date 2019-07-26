#entry problem
using TrajectoryOptimization
using Plots, LinearAlgebra



function vehicle_dynamics!(ẋ::AbstractVector, x::AbstractVector, u::Float64)
        m_cone = 200.0 #kg
        m_p = 400.0 #kg
        m = m_p + m_cone #total mass in kg
        mu = 398600.4418 #km^3.s-2
        r_G = [0.2;0.0;0.3] #COM in body frame of global spacecraft
        δ = 70*pi/180 #radians
        r_cone = 1.3 #in meters
        h_cone = r_cone/tan(δ)

        #Inertia in kg*km^2
        #J = inertia(r_cone, h_cone) #J_total at global COM
        J1 = [57.3336 0.0 10.92;
              0.0 33.3336 0.0;
              21.84 0.0 77.4]
        Jinv1 = [0.0184324 0.0 -0.00260054;
                 0.0 0.0299998 0.0;
                 -0.00520108 0.0 0.0136537]
        J = (10^(-6))*J1
        Jinv = Jinv1*(10^(-6))

        r = x[1:3]
        q = x[4:7]
        v = x[8:10]
        ω = x[11:13]

        ω_mars = [0; 0; 0.0041612*pi/180]
        v_rel = (v-cross(ω_mars, r*6378.1)) #velocity of spacecraft wrt atm

        #v_body = qrot(conj(q), v_rel)
        #α = atan(v_body[1], v_body[3])

        #control is a torque about the body z axis expressed in body frame at G global COM

        F_grav_eci = -m*mu*r/((norm(r)^3)*(6378.1^2))
        F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*6378.1, v_rel, q)
        F_aero_eci = qrot(q, F_aero_body)

        Cp = 1.0 #pitch damping coefficient
        F_total_eci = F_grav_eci + F_aero_eci
        τ_total_body = τ_aero_body - Cp*[ω[1:2]; 0.0] + [0.0;0.0;u] #computed at COM

        ẋ[1:3] = v/6378.1
        ẋ[4:7] = 0.5*qmult(q, [0; ω])
        ẋ[8:10] = F_total_eci/m
        ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
end

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
    R_e = 6378.1
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
            n̂ = [cos(v)*cos(δ); sin(v)*sin(δ); sin(δ)] #normal outward
            dA = r_i*dv*dr/sin(δ) #dA in m^2
            dA = dA*10^(-6) #dA in km^2
            r_element = [r_i*cos(v); r_i*sin(v); (r_cone-r_i)/tan(δ)]
            A = A + dA
            if n̂'*v_body > 0
                #dC = dot(n̂, v_body/norm(v_body))*dA*n̂
                #F_element = -0.5*exponential_atmosphere(h)*(norm(v_body)^2)*dC
                F_element = -0.5*exponential_atmosphere(h)*n̂'*v_body*dA*n̂*norm(v_body)*2*(n̂'*v_body/norm(v_body)) #CAREFUL
                τ_element = cross((r_element-r_G)*(10^(-3)), F_element)
                F_aero_body = F_aero_body + F_element
                τ_aero_body = τ_aero_body + τ_element
            end
        end
    end
    return F_aero_body*(A_cone/A), τ_aero_body*(A_cone/A)
end

function exponential_atmosphere(h)
    ρ0 = 1.2*10^9 #sea level density (kg/km^3)
    h0 = 8.0 #scale height (km)
    ρ = ρ0*exp(-h/h0)
end

n, m = 13, 1

entry_vehicle_model = Model(vehicle_dynamics!, n, m)

T = Float64;

#Model definition
model = entry_vehicle_model
n = model.n; m = model.m

#Initialization and goal
q0 = [1.0; 0.0; 0.0; 0.0] # unit quaternion

Re = 6378.1

x0 = zeros(T,n)
x0[1:3] = [6478.13699; 0.0; 0.0]/Re
x0[4:7] = q0
x0[8:10] = [0.0; 7.844112681678; 0.0]
x0[11:13] = [0.0; 0.0; 0.0]

xf = zero(x0)
xf[1:3] = [6478.018265205736; 39.215309132162076; -0.005906933342030818]/Re
xf[4:7] = q0
xf[8:10] = [-4.53323489922841; 3.368442794877842; -3.2449649134313976]
xf[11:13] = [0.0; 0.0; 0.0]


#Define Problem
N = 21 # number of knot points
tf = 5 # total time
dt = 0.1 # time step

U0 = [0.0 for k = 1:N-1]# initial hovering control trajectory
#Cost Function
Q = (1e-2)*Diagonal(I,n)
R = (1e-2)*Diagonal(I,m)
Qf = 1.0*Diagonal(I,n)
obj = LQRObjective(Q,R,Qf,xf,N)
#D = TrajectoryOptimization.DynamicsType

rollout!(prob)
prob = Problem(model, obj, x0=x0, integration=:rk4, N=N, dt=dt)
#initial_controls!(prob, U0); # initialize problem with controls

solve!(prob, iLQRSolverOptions{T}(verbose=true)) #solve with iLQR here
plot(prob.X)
