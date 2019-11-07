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

function sys(u, p, t)
    du = zeros(7)
    I1 = p[1] #assume principal axes
    I2 = p[2]
    I3 = p[3]
    τ = p[4:6] #see how to deal with that
    I = Diagonal([I1; I2; I3])
    I_inv = Diagonal([1/I1; 1/I2; 1/I3])
    ω = u[5:7] #angular velocity at current step
    q = u[1:4]/norm(u[1:4]) #quaternion at current step
    du[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7]))
    du[1:4] = 0.5*qmult(q, [0; ω])
    return du
end

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
x_0 = [Q; 0.01; 0.05; 0.01]
Q_0 = Matrix(Diagonal([0.001; 0.001; 0.001; 0.0001; 0.0001; 0.0001])) #covariance matrix

A1 = inv(sqrt(Q_0)); #just for computation
b1 = -A1*[0.0;0.0;0.0; x_0[5:7]];
a=1

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points = zeros(n+1, 2*n+1)
    M = -inv(A)*b #this is in dim 6
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
    for i = 1:n
        #going from dim 12 (ellipsoid) to dim 13 (propagation)
        points[1:4, 2*i-1] = qmult(C[1:4], exp_quat(V'*W[1:3, i]/2))
        #@show(exp_quat(V'*W[1:3, i]/2))
        points[5:7, 2*i-1] = C[5:7] + W[4:6, i]

        points[1:4, 2*i] = qmult(C[1:4], exp_quat(V'*(-W[1:3, i])/2))
        points[5:7, 2*i] = C[5:7] - W[4:6, i]
    end
    points[:, 2*n+1] = C
    return points #return 2n+1 points
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
    x0 = [Q; 0.01; 0.05; 0.01]
    n=6
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n+1, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
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
        x0 = x1
        #@show(X1)
    end
    return Alist, blist, centerlist, XX
end

μ = [0.0; 0.0; 0.0] #average of the perturbation (Gaussian)
Σ = [0.01 0.0 0.0; 0.0 0.01 0.0; 0.0 0.0 0.01] #covariance matrix
D = MvNormal(μ, Σ)
τ = zeros(length(μ))
rand!(D, τ)

p = [100.0; 110.0; 300.0; τ]

dt = 0.2
T = 0.0:dt:40.0
Alist, blist, centerlist, XX = propagation(A1, b1)

Plots.scatter(T, centerlist[4, :])
Plots.scatter!(T, centerlist[5, :])
Plots.scatter!(T, centerlist[6, :])
Plots.scatter!(T, centerlist[4, :])

#ref trajectory
t_span = [0.0;50.0]
dt = 0.01
y_0 = x_0
t_sim, y = rk4(sys, y_0, p, dt, t_span)


#Plot 3D ellipsoid : what should I plot ?
angles = 0.0:0.05:2*pi
angles2 = 0.0:0.05:2*pi
B = zeros(3)
b = [0.0;0.0;0.0]
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

J = 64
ellipse  = Alist[4:6, 4:6, J] \ B
plot3d!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end])

using PyPlot
gr()
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
title!("Angular velocity components")
Plots.savefig("rotation_example_angular_velocity40s_0.2step_better_integrated")
