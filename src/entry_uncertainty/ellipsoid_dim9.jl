#Trial of uncertainty propagation on translational motion AND quaternions
#only on the entry vehicle problem (angular velocity knowledge is perfect)

using Mosek
using LinearAlgebra
using Plots
using DifferentialEquations
using PyPlot
using Mosek
using MathOptInterface
using MosekTools
using JuMP
gr()

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model.jl")
include("integration.jl")
include("prop_points.jl")
include("traj_plots.jl")

a = 1

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    n = length(b)
    points2 = zeros(n, 2*n+1)
    points3 = zeros(13, 2*n+1)
    M = -inv(A)*b
    C = zeros(13)
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    C[8:10] = M[7:9]
    C[11:13] = x0[11:13]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:2*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:10, i] = points2[7:9, i]
        points3[11:13, i] = x0[11:13]
    end
    points3[:, 2*n+1] = C
    return points3
end

function points13_9(X_13, x1)
    #takes points in 13 and return points in 9
    n, m = size(X_13)
    X_9 = zeros(9, m)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    for i=1:1:m
        δq = qmult(qconj(x1[4:7]), X_13[4:7, i])
        X_9[1:3, i] = X_13[1:3, i]
        X_9[4:6, i] = 2*V*log_quat(δq)
        X_9[7:9, i] = X_13[8:10, i]
    end
    return X_9, X_13[:, end]
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

function plot_traj_center(centerlist, T)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

function prop_points_rk(X, t, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(dyna_coeffoff_COM_on_axis, X[:, i], u, 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end

function uncertainty(Alist)
    t = length(Alist[1, 1, :])
    U = zeros(t)
    for i=1:1:t
        U[i] = tr(inv(Alist[:, :, i]))
    end
    Plots.plot(U)
    Plots.xlabel!("time")
    Plots.ylabel!("uncertainty matrix trace")
end


Re = 3389.5
θ = 270*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)*1e3; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0*1e3; 0.0; 0.0; 0.0; 0.0] #okay in dimension 13: center at t=0
x0_9 = [(3389.5+125)*1e3; 0.0; 0.0; 0.0;0.0;0.0;0.0;1.0*1e3; 0.0]
Q0 = Diagonal([(100)^2;(100)^2;(100)^2;0.005^2; 0.005^2; 0.005^2;1.0;1.0;1.0]) #here we assume the change in a(vector part of the change in quaternion)

#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_9

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

δ = 70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

table_CF, table_Cτ = table_aero(δ, 0.0, 1.3, r_G)

a=1



function propagation(A1, b1)
    Δt = 100.0
    dtt = 1.0
    T = 0.0:dtt:Δt
    Re = 3389.5*1e3
    θ = 270*pi/180 #rotation angle about z-axis
    M = [-sin(θ) cos(θ) 0.0;
         0.0 0.0 1.0;
         cos(θ) sin(θ) 0.0]
    Q = mat2quat(M)
    Q = qconj(Q)
    x0 = [(3389.5+125)*1e3; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0*1e3; 0.0; 0.0; 0.0; 0.001]
    n=9
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = prop_points_rk(X1, t, dtt, u, w)
        x1 = X2[:, end]
        X_9, x0 = points13_9(X2, x1)
        A2, b2 = points2ellipse_mosek(X_9)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        x0 = x1
    end
    return Alist, blist, centerlist, XX, T
end

function X_lims(X)
    n, m = size(X)
    lim = zeros(13, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
end

function var_points(X)
    n, m = size(X)
    X_avg = (1/m)*sum(X[:, i] for i=1:1:m)
    var = (1/m)*(sum((X[:, i]-X_avg)*(X[:, i]-X_avg)' for i=1:1:m))
    sig = sqrt(var)
    return var
end

Alist, blist, centerlist, XX, T = propagation(A1, b1)

plot_traj_center(centerlist, T)
Plots.plot(centerlist[9, :])
X_lims(XX[:, :, end])
uncertainty(Alist)

var = var_points(XX[:, :, end])

oe2eci([3389.5+125; 0.000001; 45.0; 1.0; 1.0; 40.0])
x0s = [2623.5; 1676.43; 1630.38; Q; -3.32; 0.82; 0.86; 0.0; 0.0; 0.0]

x0 = [(3389.5+125)/3389.5; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; -1.0; 2.0; 0.0; 0.0; 0.0; 0.0]
t_sim4, Z4 = rk4(dyna_coeffoff_COM_on_axis, x0, [0.0], 0.01, [0.0, 150.0])

plot_ang_vel(Z4, t_sim4)
plot_altitude(Z4, t_sim4)
plot_quaternions(Z4)
plot_traj(Z4)

plot_attack_angle(Z4, t_sim4)
