using LinearAlgebra
using DifferentialEquations
using Plots
using Distributions
using Random
using Convex
using SCS
using Mosek
using JuMP
using MosekTools
using MathOptInterface
using ForwardDiff
gr()
pyplot()

include("quaternions.jl")

a=1

function sys(u, p, t)
    du = zeros(7)
    I1 = p[1] #assume principal axes
    I2 = p[2]
    I3 = p[3]
    τ = F(t) #see how to deal with that
    #τ = [0.0;0.0;0.0]
    #@show(τ)
    I = Diagonal([I1; I2; I3])
    I_inv = Diagonal([1/I1; 1/I2; 1/I3])
    ω = u[5:7] #angular velocity at current step
    q = u[1:4]/norm(u[1:4]) #quaternion at current step
    #controller
    q_ref = [0.0; 1.0; 0.0; 0.0]
    kd = 30.0 #30.0
    kp = -20.0 #20.0
    q_err = qmult(qconj(u[1:4]), q_ref) #perfect measurements
    τ_c = -kd*(u[5:7])-kp*q_err[2:4]
    #τ_c = [0.0;0.0;0.0]
    #update
    du[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7])+τ_c)
    du[1:4] = 0.5*qmult(q, [0; ω])
    return du
end
a=1
function sys!(u)
    I1 = p[1] #assume principal axes
    I2 = p[2]
    I3 = p[3]
    τ = [0.0;0.0;0.0]
    #@show(τ)
    I = Matrix(Diagonal([I1; I2; I3]))
    I_inv = Matrix(Diagonal([1/I1; 1/I2; 1/I3]))
    ω = u[5:7] #angular velocity at current step
    q = u[1:4]/norm(u[1:4]) #quaternion at current step
    #controller
    q_ref = [0.0; 1.0; 0.0; 0.0]
    kd = 30.0 #30.0
    kp = -20 #20.0
    q_err = qmult(qconj(u[1:4]), q_ref) #perfect measurements
    τ_c = -kd*u[5:7]-kp*q_err[2:4]
    #τ_c = [0.0;0.0;0.0]
    #update
    u[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7])+τ_c)
    u[1:4] = 0.5*qmult(u[1:4], [0; u[5:7]])
    return u
end

function uncertainty(Alist)
    t = length(Alist[1, 1, :])
    U = zeros(t)
    for i=1:1:t
        U[i] = tr(Alist[:, :, i])
    end
    Plots.plot(U)
    Plots.xlabel!("time")
    Plots.ylabel!("uncertainty matrix trace")
end

a=1

function rk4(f, y_0, p, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(y_star, p, t)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(y1, p, t+dt/2)
        y2 = y_star+k2*dt/2
        k3 = f(y2, p, t+dt/2)
        y3 = y_star+k3*dt
        k4 = f(y3, p, t+dt)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y
end

#=time
t_span = [0.0; 500.0]
dt = 0.05
T = t_span[1]:dt:t_span[end]

μ = [0.0; 0.0; 0.0] #average of the perturbation (Gaussian)
Σ = [0.01 0.0 0.0; 0.0 0.01 0.0; 0.0 0.0 0.01] #covariance matrix
D = MvNormal(μ, Σ)
τ = zeros(length(μ))
rand!(D, τ)

p = [100.0; 110.0; 300.0; τ]

y_0 = [1.0; 0.0; 0.0; 0.0; 0.1; 0.05; 0.01]
t_sim, y = rk4(sys, y_0, p, dt, t_span)=#


#first trial : only consider uncertaintyb in the state, will see later for τ
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x_0 = [Q; 0.2; 0.1; 0.5]
Q_0 = Matrix(Diagonal([0.01^2; 0.01^2; 0.01^2; 0.001^2; 0.001^2; 0.001^2])) #covariance matrix
V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
E0 = Matrix(Diagonal([exp_quat(V'*[0.05^2; 0.05^2; 0.05^2]/2); 0.01^2; 0.01^2; 0.01^2]))

A1 = inv(sqrt(Q_0)); #just for computation
#A1 = Matrix(Diagonal([0.01;0.01;0.01;0.05;0.05;0.05]))
b1 = -A1*[0.0;0.0;0.0; x_0[5:7]];
a=1

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points = zeros(n+1, 2*n+1)
    points2 = zeros(n, 2*n)
    points3 = zeros(n+1, 2*n+1)
    M = -inv(A)*b #this is in dim 6
    #@show(M)
    C = zeros(7, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:4] = qmult(x0[1:4], exp_quat(V'*M[1:3]/2)) #still need the reference
    #C[1:4] = [1.0; 0.0; 0.0; 0.0]
    C[5:7] = M[4:6]
    #@show(M)
    #@show(C)
    #=y = zeros(length(x)-1)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    δq = qmult(qconj(x[1:4]), x[1:4])
    if norm(δq) != 0.0
        y[1:3] = 2*V*log_quat(δq)
        y[4:6] = x[5:7]
    else
        y[1:3] = [0.0;0.0;0.0]
        y[4:6] = x[5:7]
    end
    n = length(y)
    points = zeros(n+1, 2*n+1) =#
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)

    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:2*n
        points3[1:4, i] = qmult(x0[1:4], exp_quat(V'*points2[1:3, i]/2))
        points3[5:7, i] = points2[4:6, i]
        #@show(points3[:, i])
    end
    #=
    for i = 1:n
        #going from dim 12 (ellipsoid) to dim 13 (propagation)
        points[1:4, 2*i-1] = qmult(C[1:4], exp_quat(V'*W[1:3, i]/2))
        #@show(exp_quat(V'*W[1:3, i]/2))
        points[5:7, 2*i-1] = C[5:7] + W[4:6, i]

        points[1:4, 2*i] = qmult(C[1:4], exp_quat(V'*(-W[1:3, i])/2))
        points[5:7, 2*i] = C[5:7] - W[4:6, i]
    end
    points[:, 2*n+1] = C =#
    points3[:, 2*n+1] = C
    return points3, W #return 2n+1 points
end

ellipse2points(A1, b1, x_0)
#inv(A1)

function prop_points(X, dt)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(sys, X[:, i], p, 0.001, [0.0, dt]) #because autonomous sytem
        Xnew[:, i] = Z[end, :]
    end
    return Xnew
end

function points7_6(X, x1)
    n, m = size(X);
    X_6 = zeros(n-1, m)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    for i=1:m
        δq = qmult(qconj(x1[1:4]), X[1:4, i])
        X_6[1:3, i] = 2*V*log_quat(δq)
        X_6[4:6, i] = X[5:7, i]
    end
    return X_6
end

function points2ellipse(X)
    #X is the set of points propagated through the non linear dynamics
    n, m = size(X);
    #Define Convex Optimization Problem using Convex.jl
    A = Semidefinite(n)
    b = Variable(n)
    problem = maximize(logdet(A), vcat([norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m], [A[k, j]==A[j, k] for k=1:n for j=1:n])) #[A\Matrix{Float64}(I, n, n)]))
    Convex.solve!(problem, SCSSolver(verbose = true, max_iters = 10000))
    b = b.value
    A = A.value
    @show(problem.optval)
    return A, b
end

a=1

function points2ellipse_mosek(X)
    n, m = size(X);
    s = MathOptInterface.LogDetConeTriangle(n)
    model = Model(with_optimizer(Mosek.Optimizer))
    @variable(model, A[1:n, 1:n], PSD)
    @variable(model, b[1:n])
    @variable(model , t)
    @objective(model, Max, t)
    @constraint(model, con[i = 1:m], [1.0; A*X[:, i]+b] in SecondOrderCone())
    V = [A[i, j] for j in 1:1:n for i in 1:1:j] #vectorize form of the matrix
    @constraint(model, [t;1.0;V] in s)
    #@show(con)
    #MathOptInterface.TimeLimitSec() = 0.5
    JuMP.optimize!(model)
    a = JuMP.termination_status(model)
    @show(a)
    @show(objective_value(model))
    A = JuMP.value.(A)
    b = JuMP.value.(b)
    return A, b
end

a=1

function propagation(A1, b1)
    θ = 91*pi/180 #rotation angle about z-axis
    M = [-sin(θ) cos(θ) 0.0;
         0.0 0.0 1.0;
         cos(θ) sin(θ) 0.0]
    Q = mat2quat(M)
    Q = qconj(Q)
    x0 = [Q; 0.2; 0.1; 0.5]
    n=6
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n+1, 2*n+1, length(T))
    WW = zeros(n, n, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        #@show(x0)
        X1, W = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        #@show(X1)
        X2 = prop_points(X1, dt)
        x1 = X2[:, end] #we store the last (previous center propagated)
        #@show(X2)
        X3 = points7_6(X2, x1)
        #@show(X3)
        A2, b2 = points2ellipse_mosek(X3)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        WW[:, :, i] = W
        x0 = x1
        #@show(X1)
    end
    return Alist, blist, centerlist, XX, WW
end


μ = [0.0; 0.0; 0.0] #average of the perturbation (Gaussian)
Σ = [1000.0 0.0 0.0; 0.0 1000.0 0.0; 0.0 0.0 1000.0] #covariance matrix
D = MvNormal(μ, Σ)
τ = zeros(length(μ))
rand!(D, τ)

F(t) = rand!(D, τ)

#p = [100.0; 110.0; 300.0; -0.1319; 0.08; 0.0455]
p = [100.0; 110.0; 300.0]

dt = 0.5
T = 0.0:dt:50.0
Alist, blist, centerlist, XX, WW = propagation(A1, b1)

#plot uncertainty evolution (trace)
uncertainty(Alist)

Plots.scatter(T, centerlist[1, :])
Plots.scatter!(T, centerlist[2, :])
Plots.scatter!(T, centerlist[3, :])
Plots.scatter!(T, centerlist[6, :])


#ref trajectory
t_span = [0.0;50.0]
dt = 0.01
y_0 = x_0
t_sim, y = rk4(sys, y_0, p, dt, t_span)

t_sim, y1 = rk4(sys, y_0, p, dt, t_span)

function linearize(t_sim, y)
    A = zeros(7, 7, length(t_sim))
    for i=1:1:length(t_sim)
        a = y[i, :]
        A[:, :, i] = ForwardDiff.jacobian(u ->sys!(u), a)
    end
    return A
end

A_lin = linearize(t_sim, y)

function linear_ellipse_prop(A_lin, E0)
    E = zeros(7, 7, length(t_sim))
    e = E0
    for i=1:1:length(t_sim)
        e1 = A_lin[:, :, i]*e*A_lin[:, :, i]'
        E[:, :, i] =
        e = e1
    end
    return E
end

BB = linear_ellipse_prop(A_lin, E0)

uncertainty(BB)

#Plot 3D ellipsoid
J = 60
angles = 0.0:0.05:2*pi
angles2 = 0.0:0.05:2*pi
B = zeros(3)
b = [0.0;0.0;0.0]
#b = blist[1:3, J]
function job()
      B = zeros(3)
      for i = 1:1:length(angles)
            for j = 1:1:length(angles2)
                  B = hcat(B, [cos(angles[i])*sin(angles2[j]) - b[1]; sin(angles[i])*sin(angles2[j]) - b[2]; cos(angles2[j])-b[3]])
            end
      end
      return B
end
B = job()

ellipse  = Alist[1:3, 1:3, J] \ B
Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end])
scatter3d!(WW[1, :, J], WW[2, :, J], WW[3, :, J], markersize=10.0)

J = 40
angles = 0.0:0.05:2*pi
angles2 = 0.0:0.05:2*pi
B = zeros(3)
b = [0.0;0.0;0.0]
#b = blist[4:6, J]
function job()
      B = zeros(3)
      for i = 1:1:length(angles)
            for j = 1:1:length(angles2)
                  B = hcat(B, [cos(angles[i])*sin(angles2[j]) - b[1]; sin(angles[i])*sin(angles2[j]) - b[2]; cos(angles2[j])-b[3]])
            end
      end
      return B
end
B = job()

ellipse  = Alist[4:6, 4:6, J] \ B
Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end])
scatter3d!(XX[5, :, J], XX[6, :, J], XX[7, :, J])


#Presentation plots

anim = @animate for j=1:1:100
    Plots.plot(t_sim, y[:, 5], legend = false)
    plot!(t_sim, y[:, 6])
    plot!(t_sim, y[:, 7])
    Plots.plot!(t_sim[1:50*j], y1[1:50*j, 5])
    plot!(t_sim[1:50*j], y1[1:50*j, 6])
    plot!(t_sim[1:50*j], y1[1:50*j, 7])
    Plots.scatter!(T[1:j], centerlist[4, 1:j])
    Plots.scatter!(T[1:j], centerlist[5, 1:j])
    Plots.scatter!(T[1:j], centerlist[6, 1:j])
    xlabel!("time [s]")
    ylabel!("angular velocity [rad.s-1]")
    title!("Angular velocity evolution step=$(j)")
end
gif(anim, "test_rigid_body.gif", fps = 3)


anim = @animate for J=1:1:100
    angles = 0.0:0.05:2*pi
    angles2 = 0.0:0.05:2*pi
    B = zeros(3)
    b = [0.0;0.0;0.0]
    #b = blist[4:6, J]
    function job()
          B = zeros(3)
          for i = 1:1:length(angles)
                for j = 1:1:length(angles2)
                      B = hcat(B, [cos(angles[i])*sin(angles2[j]) - b[1]; sin(angles[i])*sin(angles2[j]) - b[2]; cos(angles2[j])-b[3]])
                end
          end
          return B
    end
    B = job()
    ellipse  = Alist[4:6, 4:6, J] \ B
    if J == 1
        Plots.plot(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend = false)
    else
        Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend = false)
    end
    #scatter3d!(XX[5, :, J], XX[6, :, J], XX[7, :, J])
    title!("Ellipsoid Velocity Space Projection step=$(J)")
end
gif(anim, "test_rigid_body_ellipsoid.gif", fps = 3)

#axis-angle ellipsoids
anim = @animate for J=1:1:100
    angles = 0.0:0.05:2*pi
    angles2 = 0.0:0.05:2*pi
    B = zeros(3)
    b = [0.0;0.0;0.0]
    #b = blist[4:6, J]
    function job()
          B = zeros(3)
          for i = 1:1:length(angles)
                for j = 1:1:length(angles2)
                      B = hcat(B, [cos(angles[i])*sin(angles2[j]) - b[1]; sin(angles[i])*sin(angles2[j]) - b[2]; cos(angles2[j])-b[3]])
                end
          end
          return B
    end
    B = job()
    ellipse  = Alist[1:3, 1:3, J] \ B
    if J == 1
        Plots.plot(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend = false)
    else
        Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend = false)
    end
    #scatter3d!(XX[5, :, J], XX[6, :, J], XX[7, :, J])
    title!("Ellipsoid Attitude Space Projection step=$(J)")
end
gif(anim, "test_rigid_body_ellipsoid_attitude.gif", fps = 3)

#trace uncertainty animation

U = [tr(Alist[:, :, j]) for j=1:1:length(Alist[1, 1, :])]

anim = @animate for J=1:1:100
    Plots.plot(U[1:J], legend = false)
    Plots.xlabel!("time")
    Plots.ylabel!("uncertainty matrix trace")
end
gif(anim, "test_rigid_body_trace.gif", fps = 3)





using Plots
using PyPlot
gr()
pyplot()
using3D()
pygui(true)

fig = figure()
ax = fig[:gca](projection="3d")
Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end])

#Plot state
Plots.plot(t_sim, y[:, 1])
plot!(t_sim, y[:, 2])
plot!(t_sim, y[:, 3])
plot!(t_sim, y[:, 4])

Plots.plot(t_sim, y[:, 5])
plot!(t_sim, y[:, 6])
plot!(t_sim, y[:, 7])

#=RK4 test
function f(y, t)
    return -2*y
end

#test for RK4
t_span = [0.0; 5.0]
y_0 = 3.0
dt = 0.6
t_sim, y = rk4(f, y_0, dt, t_span)
plot(t_sim, y)
T = 0.0:0.01:5.0
plot!(T, [3.0*exp(-2*t) for t in T]) =#

xlabel!("time [s]")
ylabel!("angular vel [rad.s-1]")
title!("Angular velocity component")
Plots.savefig("rotation_example_50s_ang_vel")

#Mosek (3.5)
X1_mosek = [-0.443075 -0.524706 -0.643162 0.105911 -0.64506 0.0873887 -0.509194 -0.461816 -0.350786 -0.587458 -0.533744 -0.432474 -0.486224; 0.527601 0.62561 0.715278 -0.0753458 0.722468 -0.0582082 0.604931 0.552131 0.434497 0.683503 0.629424 0.521925 0.579388; 0.597624 0.467107 -0.251635 0.842631 -0.229813 0.843264 0.498228 0.570057 0.688599 0.343632 0.45749 0.605522 0.534933; 0.410082 0.339275 -0.106765 0.522569 -0.0955169 0.527144 0.355738 0.396122 0.462593 0.263882 0.331129 0.417018... 0.376487; -0.011734 0.0256657 0.145338 -0.131406 0.141777 -0.127845 0.0205452 -0.00661347 -0.0297794 0.0437111 0.0225832 -0.00865148 0.00696585; -0.126395 -0.17492 -0.426494 0.125179 -0.413803 0.112488 -0.164424 -0.13689 -0.0721017 -0.229213 -0.181852 -0.119462 -0.150657; 0.449384 0.524567 0.871482 0.102469 0.862256 0.111695 0.507791 0.46616 0.383007 0.590944 0.532852 0.441099 0.486976]

#
X1 = [-0.484266 -0.478715 -0.48424 -0.47874 -0.487235 -0.475652 -0.481041 -0.48195 -0.480651 -0.48234 -0.48254 -0.480448 -0.481496; 0.572648 0.577081 0.573128 0.576601 0.582606 0.567012 0.575403 0.574339 0.575094 0.574647 0.576125 0.573613 0.574871; 0.542918 0.53954 0.538774 0.543683 0.532659 0.549694 0.540772 0.541697 0.541972 0.540497 0.539907 0.54256 0.541235; 0.37789 0.383022 0.383091 0.377821 0.373425 0.387413 0.380889 0.380031 0.380142 0.380778 0.379125 0.381794 0.38046; 0.00377457... 0.00747819 0.00566794 0.00558482 0.00597458 0.00527818 0.00980094 0.00145182 0.00559188 0.00566089 0.00557441 0.00567835 0.00562638; -0.148278 -0.147716 -0.150157 -0.145837 -0.148939 -0.147055 -0.148031 -0.147962 -0.143689 -0.152305 -0.148089 -0.147905 -0.147997; 0.483509 0.482653 0.483001 0.483161 0.488046 0.478115 0.483029 0.483133 0.482988 0.483173 0.48638 0.479782 0.483081]

X2 = prop_points(X1, dt)
x1 = X2[:, end]

X2_mosek = prop_points(X1_mosek, dt)
x1_mosek = X2_mosek[:, end]

X3_mosek = points7_6(X2_mosek, x1_mosek)
X3 = points7_6(X2, x1)

(X3-X3_mosek)[:, 13]
maximum(X3-X3_mosek)

A_m, b_m = points2ellipse_mosek(X3_mosek)

A, b = points2ellipse(X3)
