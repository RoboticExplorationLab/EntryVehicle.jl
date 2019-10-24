#12 dimension ellipsoids
using Convex
using SCS
using Mosek
using LinearAlgebra
using Plots
using DifferentialEquations
using ODE
using PyPlot
pyplot()

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model_dis.jl")
include("prop_points.jl")

function ellipse2points(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    points = zeros(n+1, 2*n+1)
    M = -inv(A)*b #center of the ellipse considered
    #compute center in 13 dim
    C = zeros(length(M)+1)
    C[1:3] = M[1:3]
    #@show(M[4:6])
    C[4] = sqrt(1-M[4:6]'*M[4:6])
    C[5:end] = M[4:end]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(b)
    for i = 1:n
        #going from dim 12 (ellipsoid) to dim 13 (propagation)
        points[1:3, 2*i-1] = C[1:3] + W[1:3, i] #position
        points[4:7, 2*i-1] = qmult(C[4:7], [sqrt(abs(1-W[4:6, i]'*W[4:6, i])); W[4:6, i]]) #quat (eigenvalues ?)
        points[8:13, 2*i-1] = C[8:13] + W[7:12, i] #remaining

        points[1:3, 2*i] = C[1:3] - W[1:3, i] #position
        points[4:7, 2*i] = qmult(C[4:7], [sqrt(abs(1-W[4:6, i]'*W[4:6, i])); -W[4:6, i]]) #quat (eigenvalues ?)
        points[8:13, 2*i] = C[8:13] - W[7:12, i] #remaining
    end
    points[:, 2*n+1] = C
    return points #return 2n+1 points (in dim 13 from ellipsoid in 12) so far
end

function points2ellipse(X)
    #X is the set of points propagated through the non linear dynamics
    n, m = size(X);
    #Define Convex Optimization Problem using Convex.jl
    A = Semidefinite(n)
    b = Variable(n)
    problem = maximize(logdet(A), vcat([norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m], [A[k, j]==A[j, k] for k=1:n for j=1:n])) #[A\Matrix{Float64}(I, n, n)]))
    Convex.solve!(problem, SCSSolver(verbose = false))
    b = b.value
    A = A.value
    return A, b
end

function points13_12(X)
    n, m = size(X);
    X_12 = zeros(n-1, m)
    for i=1:m
        X_12[1:3, i] = X[1:3, i]
        X_12[4:6, i] = X[5:7, i]
        X_12[7:12, i] = X[8:13, i]
    end
    return X_12
end

function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

#should start with ellipse in dim 12 then I believe
#=function prop
#X1 = ellipse2points(A1, b1) #input in dim 12 and output in dim 13
#X2 = propagation(X1) #input in dim 13 and output in dim 13
#X3 = points13_12(X2) #input in dim 13 and output in dim 12
#A2, b2 = points2ellipse(X3) #ouput ellipse in dim 12 =#


Re = 3389.5
θ = 110*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal(0.000000001*ones(12))
#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A0 = inv(sqrt(Q0))
b0 = -A1*x0

Δt = 3.0 #length simulation
dt = 0.1
T = 0.0:dt:Δt

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

n = length(x0)
δ = 70*pi/180
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_G)

function propagation(A1, b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n+1, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1) #return a set of points
        X2 = prop_points_continuous(X1, dt, u, w) #prop_points(X1, t, dt, u, w)
        X3 = points13_12(X2)
        A2, b2 = points2ellipse(X3)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        #@show(X1)
    end
    return Alist, blist, centerlist, XX
end

Alist, blist, centerlist, XX = propagation(A0, b0)

#test plots ellipses (in X and Y) : not good
j = 31
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
Plots.plot!(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
Plots.scatter(XX[1, :, j], XX[2, :, j])

plot_traj_center(centerlist)
Plots.savefig("uncertainty-0.001-10sec-2")
xlabel!("X")
ylabel!("Y")
