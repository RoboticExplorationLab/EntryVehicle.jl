using LinearAlgebra
using DifferentialEquations
using Plots
using Distributions
using Random
using Convex
using SCS
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
    #@show(τ)
    I = Diagonal([I1; I2; I3])
    I_inv = Diagonal([1/I1; 1/I2; 1/I3])
    ω = u[5:7] #angular velocity at current step
    q = u[1:4]/norm(u[1:4]) #quaternion at current step
    #controller
    q_ref = [0.0; 1.0; 0.0; 0.0]
    kd = 30.0
    kp = -20.0
    q_err = qmult(qconj(u[1:4]), q_ref) #perfect measurements
    τ_c = -kd*u[5:7]-kp*q_err[2:4]
    #τ_c = [0.0;0.0;0.0]
    #update
    du[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7])+τ_c)
    du[1:4] = 0.5*qmult(q, [0; ω])
    return du
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
Q_0 = Matrix(Diagonal([0.01; 0.01; 0.01; 0.001; 0.001; 0.001])) #covariance matrix

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
    return A, b
end

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
        X1, W = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = prop_points(X1, dt)
        x1 = X2[:, end] #we store the last (previous center propagated)
        X3 = points7_6(X2, x1)
        A2, b2 = points2ellipse(X3)
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

Plots.scatter(T, centerlist[4, :])
Plots.scatter!(T, centerlist[5, :])
Plots.scatter!(T, centerlist[6, :])
Plots.scatter!(T, centerlist[6, :])

norm(centerlist[1:3, 82])

#ref trajectory
t_span = [0.0;50.0]
dt = 0.01
y_0 = x_0
t_sim, y = rk4(sys, y_0, p, dt, t_span)

norm(y[9000, 2:4])
norm(y[8000, 2:4])

norm(y[10001, 2:4])
inv(Alist[1:6, 1:6, 180])

Alist[1:6, 1:6, 200]
inv(Alist[1:6, 1:6, 40]*Alist[1:6, 1:6, 40])

E = [exp_quat(y[i, 1:4]) for i in 1:1:length(t_sim)]

#Plot 3D ellipsoid
J = 101
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
scatter3d!(WW[1, :, J], WW[2, :, J], WW[3, :, J], markersize=20.0)

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
Plots.plot(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end])
scatter3d!(XX[5, :, J], XX[6, :, J], XX[7, :, J])







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
Plots.savefig("rotation_example_quat_true_ellip_150s_notorque")
