#Trial of uncertainty propagation on translational motion only on the entry vehicle problem
using Mosek
using LinearAlgebra
using Plots
using DifferentialEquations
using ODE
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

a=1


function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    points2 = zeros(n, 2*n+1)
    points3 = zeros(13, 2*n+1)
    M = -inv(A)*b
    C = zeros(13)
    C[1:3] = M[1:3]
    C[4:7] = x0[4:7]
    C[8:10] = M[4:6]
    C[11:13] = x0[11:13]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:2*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = x0[4:7]
        points3[8:10, i] = points2[4:6, i]
        points3[11:13, i] = x0[11:13]
    end
    points3[:, 2*n+1] = C
    return points3
end

function points13_6(X_13)
    #takes points in 13 and return points in 6
    n, m = size(X_13)
    X_6 = zeros(6, m)
    for i=1:1:m
        X_6[1:3, i] = X_13[1:3, i]
        X_6[4:6, i] = X_13[8:10, i]
    end
    return X_6, X_13[:, 1]
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

function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

#=
function rk4(f, y_0, u, dt, t_span)
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
        k1 = f(t, y_star, u)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, u)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, u)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, u)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y
end =#

function prop_points_rk(X, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(dyna_coeffoff, X[:, i], u, 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
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
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; -2.0; 3.0; 0.0; 0.0; 0.0; 0.0] #okay in dimension 13: center at t=0
x0_6 = [(3389.5+125)/Re; 0.0; 0.0; 0.0; 1.0; 0.0]
Q0 = Diagonal(0.00000001*ones(12)) #here we assume the change in a(vector part of the change in quaternion)
Q0 = Diagonal([(0.1/Re)^2;(0.1/Re)^2;(0.1/Re)^2;(1e-5)^2;(1e-5)^2;(1e-5)^2]) #here we assume the change in a(vector part of the change in quaternion)

#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_6

Δt = 300.0 #length simulation
dt = 1.0
T = 0.0:dt:Δt
#T = t_sim_nom

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

δ = 70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

a=1

function propagation(A1, b1)
    Re = 3389.5
    θ = 91*pi/180 #rotation angle about z-axis
    M = [-sin(θ) cos(θ) 0.0;
         0.0 0.0 1.0;
         cos(θ) sin(θ) 0.0]
    Q = mat2quat(M)
    Q = qconj(Q)
    x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
    n=6
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = prop_points_rk(X1, dt, u, w)
        X_6, x0 = points13_6(X2)
        A2, b2 = points2ellipse_mosek(X_6)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX
end

function X_lims(X)
    n, m = size(X)
    lim = zeros(13, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
end

Alist, blist, centerlist, XX = propagation(A1, b1)

plot_traj_center(centerlist)
uncertainty(Alist)
lim = X_lims(XX[:, :, end])

Plots.plot(centerlist[1, :], centerlist[2,:])
Plots.plot(Z4'[1, :], Z4'[2, :])
t_sim4, Z4 = rk4(dyna_coeffoff, x0, [0.0], 0.001, [0.0, Δt])

plot_ang_vel(Z4', t_sim4)
plot_altitude(Z4', t_sim4)
plot_quaternions(Z4')
plot_traj(Z4')

plot_attack_angle(Z4', t_sim4)
plot_vel(Z4', t_sim4)

j = 1000
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]); sin(angles[i])]
end

ellipse  = Alist[1:2, 1:2, j] \ B
Plots.plot!(ellipse[1, :].+centerlist[1, j], ellipse[2, :].+centerlist[2, j])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!(XX[1, :, j], XX[2, :, j])
