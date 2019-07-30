#test and learning on Trajopt module
using TrajectoryOptimization
using Plots, LinearAlgebra


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

function exponential_atmosphere(h)
    ρ0 = 0.026455*10^9 #sea level density (kg/km^3)
    h0 = 10.8 #scale height on Mars(km)
    ρ = ρ0*exp(-h/h0)
end

#dynamics

function dyna!(ẋ, x, u)
    # ẋ = zeros(13)

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
    #=J = [57.3336 0.0 21.84;
            0.0 33.3336 0.0;
            21.84 0.0 77.4]
    Jinv = [0.0184324 0.0 -0.00520108;
             0.0 0.0299998 0.0;
             -0.00520108 0.0 0.0136537] =#

    #No mass offset
    J = [77.208 0.0 0.0;
          0.0 77.208 0.0;
          0.0 0.0 101.4]

    Jinv = [0.012952 0.0 0.0;
            0.0 0.012952 0.0;
            0.0 0.0  0.00986193]

    r = x[1:3]
    q = x[4:7]/norm(x[4:7])
    v = x[8:10]
    ω = x[11:13]

    if norm(r) > 1.0

    ω_mars = [0; 0; 0.0040612*pi/180]
    v_rel = (v) #-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm

    #v_body = qrot(conj(q), v_rel)
    #α = atan(v_body[1], v_body[3])

    #control is a torque about the body z axis expressed in body frame at G global COM

    F_grav_eci = -m*mu*r/((norm(r)^3)*(Re^2))
    F_aero_body, τ_aero_body = compute_aero(δ, r_cone, r_G, r*Re, v_rel, q)
    F_aero_eci = qrot(q, F_aero_body)

    Cp = 1.0 #pitch damping coefficient
    F_total_eci = F_grav_eci + F_aero_eci
    τ_total_body = τ_aero_body #- Cp*[ω[1:2]; 0.0] + [0.0;0.0;u[1]] #computed at COM #[0.0;0.0;0.0]

    ẋ[1:3] = v/Re
    ẋ[4:7] = 0.5*qmult(q, [0; ω])
    ẋ[8:10] = F_total_eci/m
    ẋ[11:13] = Jinv*(τ_total_body-cross(ω, J*ω))
else
    # ẋ = zeros(13)
end
end

T = Float64


n, m = 13, 1

entry_vehicle_model = Model(dyna!, n, m)

#Model definition
model = entry_vehicle_model
n = model.n; m = model.m

Re = 3389.5
#test dynamic model
q00 = [1.0; 0.0; 0.0; 0.0]
qff = [1.0; 0.0; 0.0; 0.0]
x0 = zeros(T, n)
x0[1:3] = [3389.5+125; 0.0; 0.0]/Re
x0[4:7] = q00
x0[8:10] = [0.0; 1.0; 0.0]
x0[11:13] = [0.0; 0.0; 0.0]

#expected values after 100 seconds and no control inputs
xf = zero(x0)
#xf[1:3] = [1.0082492161173011; 0.1223730525679208;-0.0003671264111757556]
xf[4:7] = q00
#xf[8:10] = [-0.9472041894063383; 7.7470514173024485; -0.04666470034355125]
xf[11:13] = [0.0; 0.0; 0.0]
#here for 550 seconds without control
#xf[1:3] = [0.8015097123420436; 0.5969821372031944; -0.034304797979191365]
#xf[8:10] = [-4.5333007135526495; 3.368324196015978; -3.2450298665674895]
#here with shift, target is before 1 degree
#xf[1:3] = [0.811231302755334; 0.5848325190529827; -0.0334603593626661]
#xf[8:10] = [-4.443867255757695;3.4421107961654713; -3.2409461334862266]
xf[1:3] = [0.9962;0.085; 0.0]
xf[8:10] = [0.0; 0.0; 0.0]

#xdot = zeros(n)
#x, u = x0, [0.0]
#evaluate!(xdot,model,x,u)

N = 101
Q = Diagonal(0.01I,n)
R = Diagonal(0.01I,m)
Qf = Diagonal(10.0I,n) #Put float everywhere

t0 = 0
tf = 500  #total time ? yeah pretty sure
obj = LQRObjective(Q, R, Qf, xf, N)


#constraint
bnd = BoundConstraint(n, m, u_min = -0.001, u_max = 0.001)
constraints = [bnd] #only bounds constraints on my system, put the []
bnd isa ConstraintSet
constraints isa ConstraintSet
CON = Constraints(constraints, N)

#=c1(v,x,u) = v[1] = u[1] - 0.001
p1 = 1
con1 = Constraint{Inequality}(c1, n, m, p1, :mycon1)
c2(v,x,u) = v[1] = -u[1] - 0.001
p2 = 2
con2 = Constraint{Inequality}(c2, n, m, p2, :mycon2)
constraints = con1 + con2
CON = Constraints(constraints, N) =#

#solver options
max_con_viol = 1.0e-10 #been changed
verbose=true


opts_ilqr = iLQRSolverOptions{T}(verbose=verbose,live_plotting=:off)

opts_al = AugmentedLagrangianSolverOptions{T}(verbose=verbose,
    opts_uncon=opts_ilqr,
    cost_tolerance=1.0e-6,
    cost_tolerance_intermediate=1.0e-2,
    constraint_tolerance=max_con_viol,
    penalty_scaling=50.,
    penalty_initial=10.)

opts_pn = ProjectedNewtonSolverOptions{T}(verbose=verbose,
    feasibility_tolerance=max_con_viol,
    solve_type=:feasible)

opts_altro = ALTROSolverOptions{T}(verbose=verbose,
    opts_al=opts_al,
    R_inf=1.0e-1,
    resolve_feasible_problem=false,
    opts_pn=opts_pn,
    projected_newton=true,
    projected_newton_tolerance=1.0e-3);

prob = Problem(model, obj, x0 = x0, integration=:rk4, constraints = CON, N=N, tf=tf)
TrajectoryOptimization.solve!(prob, opts_al)


#Plots states and control
x = zeros(N)
y = zeros(N)
z = zeros(N)
q0 = zeros(N)
q1 = zeros(N)
q2 = zeros(N)
q3 = zeros(N)
vx = zeros(N)
vy = zeros(N)
vz = zeros(N)
ωx = zeros(N)
ωy = zeros(N)
ωz = zeros(N)
alt = zeros(N)

for i in 1:1:length(prob.X)
    array = prob.X[i]
    x[i] = array[1]
    y[i] = array[2]
    z[i] = array[3]
    q0[i] = array[4]
    q1[i] = array[5]
    q2[i] = array[6]
    q3[i] = array[7]
    vx[i] = array[8]
    vy[i] = array[9]
    vz[i] = array[10]
    ωx[i] = array[11]
    ωy[i] = array[12]
    ωz[i] = array[13]
    alt[i] = sqrt(x[i]^2+y[i]^2+z[i]^2)*6378.1
end

#=
plot(q0)
plot!(q1)
plot!(q2)
plot!(q3)

plot(ωx)
plot!(ωy)
plot!(ωz)

plot(alt)
plot(prob.U) =#
