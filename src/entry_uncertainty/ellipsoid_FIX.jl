#12 dimension ellipsoids
using Convex
using SCS
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
using Distributions
using Random
gr()

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model.jl")
include("integration.jl")
include("prop_points.jl")
include("traj_plots.jl")

#used to compute the time sequence for the nominal trajectory

 ####################################
 ######Simulation Parameters#########
 ####################################

Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
      0.0 0.0 1.0;
      cos(θ) sin(θ) 0.0]
Q = mat2quat(M) #CHANGE THAT
Q = qconj(Q)
#x0 = [(3389.5+125)/Re, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Δt = 300 #length simulation

 ####################################
 #####Dynamics - Integration#########
 ####################################


 #new part for offline aerodynamic coefficients computation
δ = 70*pi/180
r_cone = 1.3
r_min = 0.2
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_min, r_G)

t_sim_nom, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)
t_sim4, Z4 = rk4(dyna_coeffoff, x0, [0.0], 0.001, [0.0, Δt])

plot_attack_angle(Z4, t_sim4)

a = 1

Plots.scatter(centerlist[11, :])
Plots.scatter(Z4[13, :])
Plots.scatter!(Z[1, :]*Re, Z[2, :]*Re)
plot_ang_vel(Z4, t_sim4)
plot_altitude(Z4, t_sim4)
plot_quaternions(Z4)
Plots.scatter!(centerlist[10, :])

######################################
########Ellipsoid propagation#########
######################################

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    n = length(b)
    points = zeros(n+1, 2*n+1)
    #M = -inv(A)*b #center of the ellipse considered in dimension 12 here not necessary here
    #compute center in 13 dim
    M = -inv(A)*b #this is in dim 12
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], [sqrt(1-M[4:6]'*M[4:6]); M[4:6]]) #still need the reference
    C[8:13] = M[7:12]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i = 1:n
        #going from dim 12 (ellipsoid) to dim 13 (propagation)
        points[1:3, 2*i-1] = C[1:3] + W[1:3, i] #position
        points[4:7, 2*i-1] = qmult(C[4:7], [sqrt(1-W[4:6, i]'*W[4:6, i]); W[4:6, i]]) #no conjugate here for sure
        points[8:13, 2*i-1] = C[8:13] + W[7:12, i] #remaining

        points[1:3, 2*i] = C[1:3] - W[1:3, i] #position
        points[4:7, 2*i] = qmult(C[4:7], [sqrt(1-W[4:6, i]'*W[4:6, i]); -W[4:6, i]])
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
    Convex.solve!(problem, SCSSolver(verbose = true, max_iters = 15000))
    b = b.value
    A = A.value
    return A, b
end

function points13_12(X, x1)
    n, m = size(X);
    X_12 = zeros(n-1, m)
    for i=1:m
        X_12[1:3, i] = X[1:3, i]
        #X_12[4:6, i] = X[5:7, i] #nah need deltaq with respct to a ref quat
        a = qmult(qconj(x1[4:7]), X[4:7, i])
        X_12[4:6] = a[2:end] #take vector part of the error. Here qconj for sure too
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
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0] #okay in dimension 13: center at t=0
x0_12 = [(3389.5+125)/Re; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0] #because the center with respect to himself is 0
Q0 = Diagonal([(0.1/Re)^2;(0.1/Re)^2; (0.1/Re)^2; 0.001^2; 0.001^2; 0.001^2; 0.05^2; 0.05^2;0.05^2; 0.0000001; 0.0000001; 0.0000001]) #here we assume the change in a(vector part of the change in quaternion)
#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_12

Δt = 30.0 #length simulation
dt = 0.1
T = 0.0:dt:Δt
#T = t_sim_nom

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

n = length(x0)
δ = 70*pi/180
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_G)
n = 12

function propagation(A1, b1)
    Re = 3389.5
    θ = 91*pi/180 #rotation angle about z-axis
    M = [-sin(θ) cos(θ) 0.0;
         0.0 0.0 1.0;
         cos(θ) sin(θ) 0.0]
    Q = mat2quat(M)
    Q = qconj(Q)
    x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
    n=12
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n+1, 2*n+1, length(T))
    for i=1:1:length(T)
        #=if i ==1
            dt = T[i]
            t = T[i]
        else
            t = T[i] #don't really need that, autonomous system here in fact
            dt = T[i]-T[i-1]
        end=#
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = prop_points_last(X1, dt, u, w)
        x1 = X2[:, end] #we store the last (previous center propagated)
        X3 = points13_12(X2, x1)
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

Alist, blist, centerlist, XX = propagation(A1, b1)

#test plots ellipses (in X and Y) : not good
j = 500
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i])-blist[1, j]; sin(angles[i])-blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
Plots.plot!(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!(XX[1, :, j], XX[2, :, j])

T = 1:100 #100 means I go for 50 ec
plot_traj_center(centerlist[:, 1:100])

Plots.savefig("Propagation-entry-50sec-0.5 step- 15000 maxiters-more uncertainty")

#Alist[:, :, 50] #check symmetry















#Test with uncertainty only on position and velocity to see if ellipsoids computation works

Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
      0.0 0.0 1.0;
      cos(θ) sin(θ) 0.0]
Q = mat2quat(M) #CHANGE THAT
Q = qconj(Q)
#x0 = [(3389.5+125)/Re, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Δt = 100 #length simulation

 ####################################
 #####Dynamics - Integration#########
 ####################################

t_sim_nom, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)

 #new part for offline aerodynamic coefficients computation
δ = 70*pi/180
r_min = 0.2
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero_spherecone(δ, r_min, r_cone, r_G)

t_sim_nom, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)

plot_traj(Z)
plot_altitude(t_sim_nom, Z)
plot_ang_vel(Z)
Plots.savefig("100s_ellipse_entry")

a=1

######################################
########Ellipsoid propagation#########
######################################

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points2 = zeros(n, 2*n)
    points3 = zeros(n+1, 2*n+1)
    M = -inv(A)*b
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    #C[1:4] = [1.0; 0.0; 0.0; 0.0]
    C[8:13] = M[7:12]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:2*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:13, i] = points2[7:12, i]
    end
    points3[:, 2*n+1] = C
    return points3, W #return 2n+1 points
end

function ellipse2points2(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points2 = zeros(n, 6*n)
    points3 = zeros(n+1, 6*n+1)
    M = -inv(A)*b
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    #C[1:4] = [1.0; 0.0; 0.0; 0.0]
    C[8:13] = M[7:12]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)
    for i =1:n
        points2[:, 6*i-5] = M + 0.5*W[:, i]
        points2[:, 6*i-4] = M - 0.5*W[:, i]
        points2[:, 6*i-3] = M - 0.8*W[:, i]
        points2[:, 6*i-2] = M + 0.8*W[:, i]
        points2[:, 6*i-1] = M - W[:, i]
        points2[:, 6*i] = M + W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:6*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:13, i] = points2[7:12, i]
    end
    points3[:, 6*n+1] = C
    return points3, W #return 2n+1 points
end

ellipse2points2(A1, b1, x_0)
ellipse2points(A1, b1, x_0)

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

function points13_12(X, x1)
    n, m = size(X);
    X_12 = zeros(n-1, m)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    for i=1:m
        δq = qmult(qconj(x1[4:7]), X[4:7, i])
        X_12[4:6, i] = 2*V*log_quat(δq)
        X_12[1:3, i] = X[1:3, i]
        X_12[7:12, i] = X[8:13, i]
    end
    return X_12
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
end

function prop_points_last(X, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(dyna_coeffoff, X[:, i], u, 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[end, :]
    end
    return Xnew
end


#should start with ellipse in dim 12 then I believe
#=function prop
#X1 = ellipse2points(A1, b1) #input in dim 12 and output in dim 13
#X2 = propagation(X1) #input in dim 13 and output in dim 13
#X3 = points13_12(X2) #input in dim 13 and output in dim 12
#A2, b2 = points2ellipse(X3) #ouput ellipse in dim 12 =#

Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0] #okay in dimension 13: center at t=0
x0_12 = [(3389.5+125)/Re; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]#because the center with respect to himself is 0
Q0 = Diagonal(0.00000001*ones(12)) #here we assume the change in a(vector part of the change in quaternion)
Q0 = Diagonal([(0.05/Re)^2;(0.05/Re)^2; (0.05/Re)^2; 0.01^2; 0.01^2; 0.01^2; (1e-5)^2; (1e-5)^2;(1e-5)^2; (1e-20); (1e-20); (1e-20)]) #here we assume the change in a(vector part of the change in quaternion)

#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_12

Δt = 200.0 #length simulation
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
n = 12

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
    n=12
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
        @show(X1)
        X2 = prop_points_last(X1, dt, u, w)
        x1 = X2[:, end] #we store the last (previous center propagated)
        @show(X2)
        X3 = points13_12(X2, x1)
        @show(X3)
        A2, b2 = points2ellipse_mosek(X3)
        @show(A2, b2)
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

Alist, blist, centerlist, XX, WW = propagation(A1, b1)
X1[:, 22]


function X_lims(X)
    n, m = size(X)
    lim = zeros(13, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
end

t_sim, Z = rk4(dyna_coeffoff, X1[:, 22], u, 0.01, [0.0, 300])
plot_traj(Z')
plot_altitude(Z', t_sim)

#q1 = [-0.0207074, -0.894882, -0.04728, 0.443307]
#norm(q1)

X1 = [1.01222 1.01132 1.01151 1.01203 1.01151 1.01203 1.01154 1.012 1.01098 1.01257 1.01323 1.01031 1.01199 1.01156 1.00916 1.01439 1.00233 1.02121 1.01173 1.01181 0.696013 1.32753 1.01169 1.01185 1.01177; 0.0638269 0.0643482 0.0645811 0.063594 0.0646373 0.0635378 0.0646608 0.0635143 0.0650376 0.0631375 0.0626791 0.065496 0.0638332 0.0643419 0.0668646 0.0613105 0.0744262 0.0537489 0.0641443 0.0640308 0.372615 -0.24444 0.0642305 0.0639446 0.0640876; -0.00108567 -0.000563123 -0.00027466 0.00137413 0.00035366 -0.00200246 0.000312324 -0.00196112 0.000708098 -0.00235689 -0.00238755 0.000738758 -0.00102544 -0.000623357 0.00253669 -0.00418548 0.0121953 -0.0138441 -0.00073338 -0.000915415 0.332424 -0.334073 -0.000582535 -0.00106626 -0.000824397; -0.762057 -0.762756 -0.762846 -0.761967 -0.7631 -0.761713 0.714412 0.763704 0.873331 0.498849 0.809106 0.749055 -0.76226 -0.762553 -0.768479 -0.756167 -0.780622 -0.742347 -0.75847 -0.766274 0.0421331 0.138771 -0.763116 -0.761684 -0.762407; 0.15789 0.157124 0.157037 0.157977 0.156874 0.15814 -0.339651 -0.236324 -0.316784 -0.148453 0.00461251 -0.270045 0.157312 0.157701 0.152355 0.162617 0.136781 0.177809 0.161576 0.153412 -0.313286 0.194683 0.153775 0.161234 0.157507; 0.615702 0.61528 0.615173 0.615809 0.614897 0.616084 -0.554768 -0.509139 -0.300136 -0.787239 -0.586923 -0.526384 0.615995 0.614985 0.611366 0.619505 0.603784 0.625971 0.620471 0.610454 0.151644 -0.233439 0.615176 0.615795 0.615491; -0.123498 -0.122253 -0.122346 0.123405 -0.122359 -0.123391 0.257848 0.318881 0.21647 0.330707 -0.0291032 0.298181 -0.121505 -0.124246 -0.111623 -0.134095 -0.0858229 -0.159576 -0.116759 -0.128973 0.936526 -0.942522 -0.12476 -0.120986 -0.122875; -0.768587 -0.769016 -0.769056 -0.768547 -0.769003 -0.768601 -0.770352 -0.767251 -0.769105 -0.768499 -0.771311 -0.766292 -0.763186 -0.774418 -0.767472 -0.770131 -0.760665 -0.776939 -0.768678 -0.768925 -0.320008 -1.2176 -0.769371 -0.768232 -0.768802; 0.974649 0.979877 0.98004 0.974486 0.980624 0.973902 0.978795 0.975731 0.992912 0.961613 0.95301 1.00152 0.978592 0.975933 1.02033 0.934195 1.12373 0.830795 0.978661 0.975865 6.05488 -4.10035 0.976705 0.97782 0.977263; 0.00591342 0.0247978 0.0256942 0.00501696 0.0283753 0.00233592 0.0240014 0.00670973 0.0599629 -0.0292517 -0.0688514 0.0995626 0.0234926 0.00721853 0.161824 -0.131113 0.537783 -0.507072 0.0191541 0.011557 17.8712 -17.8405 0.0150581 0.0156531 0.0153556; 0.00154332 0.00162347 0.00164015 0.00152664 0.00167441 0.00149238 -0.0162106 0.0193774 -0.00442486 0.00759165 -0.00383979 0.00700658 0.00170666 0.00146013 0.00298122 0.000185573 0.00538196 -0.00221517 0.417064 -0.413897 0.00235066 0.000816129 -0.127178 0.130345 0.00158339; 6.88498 7.5165 7.50927 6.89221 7.53399 6.86749 6.00967 8.39182 7.38386 7.01763 4.07269 10.3288 7.64953 6.75195 12.2784 2.12312 25.0566 -10.6551 7.20151 7.19997 860.857 -846.455 7.19971 7.20177 7.20074; -3.7289e-5 0.000118321 0.000183488 -0.000102457 0.000282378 -0.000201347 0.0082584 -0.00817737 -0.00198248 0.00206351 -0.00134146 0.00142249 -0.000528906 0.000609938 -0.000516872 0.000597904 -0.000256985 0.000338016 -0.128721 0.128802 -0.000987342 0.00106837 0.208836 -0.208755 4.05158e-5]


X1 = [1.20177 0.890677 1.07675 1.0157 1.04705 1.0454 1.02908 1.06336 1.08187 1.01058 1.09908 0.993362 4.49509 -2.40265 1.89903 0.193416 1.06158 1.03086 1.04622 1.04623 -1.83336 3.92581 1.04623 1.04622 1.04622; 0.0783621 0.0173119 0.0562627 0.0394114 0.0479876 0.0476864 0.0419948 0.0536793 0.053402 0.0422721 0.0610439 0.0346302 0.698045 -0.602371 0.225937 -0.130263 0.0499044 0.0457697 0.0478569 0.0478171 -0.532628 0.628302 0.0478279 0.0478462 0.047837; 0.000923017 -0.000728721 0.000247774 -5.34784e-5 0.000168777 2.55187e-5 4.04772e-5 0.000153819 0.000302053 -0.000107757 0.000336584 -0.000142288 0.0189147 -0.0187204 0.00454049 -0.0043462 0.000222803 -2.85074e-5 9.18149e-5 0.000102481 -0.0153107 0.015505 9.86378e-5 9.5658e-5 9.71479e-5; 0.126955 0.133995 0.12841 0.13264 0.130528 0.130527 -0.116198 0.0977282 -0.0687724 -0.230144 -0.161794 0.17031 0.347312 -0.136176 -0.297825 0.463546 0.137353 0.12369 0.130527 0.130528 0.199529 -0.114963 0.130531 0.130524 0.130528; -0.198993 -0.136056 -0.174523 -0.16066 -0.167749 -0.167439 0.0431416 -0.101409 0.0110405 0.203031 -0.120028 0.439708 -0.563772 0.286212 -0.798313 0.525177 -0.161924 -0.173259 -0.167602 -0.167586 -0.875424 0.614149 -0.167586 -0.167602 -0.167594; -0.516171 -0.513774 -0.516571 -0.514081 -0.515339 -0.515354 0.882663 0.120414 0.769667 0.766283 0.789136 0.650006 -0.248356 -0.589502 -0.410263 -0.0220965 -0.512281 -0.518348 -0.515344 -0.515349 -0.326544 0.522184 -0.515346 -0.515347 -0.515346; 0.823246 0.836333 0.828305 0.831985 0.830139 0.830192 -0.453238 -0.982623 -0.634542 -0.564363 -0.58014 -0.595841 0.706921 0.742906 0.324906 0.713237 0.832085 0.828179 0.830166 0.830166 -0.295081 0.580354 0.830167 0.830164 0.830166; 3.20801 -3.68974 0.409341 -0.891076 -0.22205 -0.259685 -0.0760178 -0.405717 0.781004 -1.26274 0.149228 -0.630963 82.2702 -82.7519 16.5008 -16.9826 0.171149 -0.652884 -0.240858 -0.240877 -62.2099 61.7282 -0.240873 -0.240862 -0.240868; 1.9617 0.256086 1.28699 0.930792 1.11333 1.10445 0.278514 1.93927 1.38738 0.830404 2.85828 -0.640494 17.8506 -15.6328 10.8145 -8.59675 1.12857 1.08921 1.10888 1.1089 -23.0106 25.2284 1.1089 1.10889 1.10889; 0.0129611 -0.017764 -0.000334048 -0.00446878 -0.00227576 -0.00252707 0.00621259 -0.0110154 0.00145716 -0.00626 -0.0195202 0.0147173 0.409615 -0.414418 0.0172816 -0.0220844 0.0173115 -0.0221143 -0.00239759 -0.00240524 -0.224769 0.219966 -0.00241468 -0.00238816 -0.00240142; -3.86459e-6 6.18738e-6 2.10733e-5 -1.87505e-5 -4.17159e-6 6.49438e-6 2.82506e-6 -5.02271e-7 1.13317e-5 -9.00887e-6 1.41543e-5 -1.18315e-5 1.01838e-5 -7.86101e-6 -9.97071e-6 1.22935e-5 4.98401e-6 -2.66122e-6 0.000669072 -0.00066675 1.11241e-6 1.21038e-6 -0.000271304 0.000273627 1.1614e-6; 42.5381 48.2973 44.8373 45.9982 45.4023 45.4331 43.3251 47.5103 47.1698 43.6656 46.7081 44.1273 -16.5513 107.387 21.2982 69.5372 45.1954 45.6401 45.4177 45.4177 896.65 -805.814 45.4177 45.4177 45.4177; 1.8802e-6 -7.87899e-6 -1.21321e-5 6.13327e-6 -1.50951e-6 -4.48928e-6 1.87736e-6 -7.87616e-6 -7.73593e-6 1.73713e-6 -1.94324e-5 1.34336e-5 -8.05561e-6 2.05681e-6 3.16571e-6 -9.16451e-6 -1.62589e-5 1.02601e-5 -0.000275465 0.000269466 -2.88156e-6 -3.11724e-6 0.00043558 -0.000441578 -2.9994e-6]

X1[:, 1]
X1[:, 23]
X1[:, 21]
X1[:, 1]



function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]*Re
        Y[j] = centerlist[2, j]*Re
    end
    Plots.scatter(X, Y)
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

plot_traj_center(centerlist)
uncertainty(Alist)

Plots.scatter(centerlist[6, :])

maximum(centerlist[1, :])
minimum(centerlist[1, :])

a=1

#Plot 3D ellipsoid
anim = @animate for J=1:5:100
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

    ellipse  = Alist[10:12, 10:12, J] \ B
    if J ==1
        Plots.plot(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend =false)
    else
        Plots.plot!(ellipse[1, 2:end], ellipse[2, 2:end], ellipse[3, 2:end], legend =false)
    end
    #xlims!((-0.03, 0.03))
    ylims!((-0.0004, 0.0004))
    #zlims!((-0.02, 0.02))
    title!("Ellipsoid ang vel step=$(J)")
end
gif(anim, "plot_ellipse_ang_vel.gif", fps = 4)

scatter3d!(WW[1, :, J], WW[2, :, J], WW[3, :, J], markersize=10.0)


#MC on entry vehicle
Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0] #okay in dimension 13: center at t=0
x0_12 = [(3389.5+125)/Re; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]#because the center with respect to himself is 0


D1 = Uniform(x_0_6[1]-σ_x, x_0_6[1]+σ_x)
D2 = Uniform(x_0_6[2]-σ_y, x_0_6[2]+σ_y)
D3 = Uniform(x_0_6[3]-σ_z, x_0_6[3]+σ_z)
D4 = Uniform(x_0_6[4]-σ_ωx, x_0_6[4]+σ_ωx)
D5 = Uniform(x_0_6[5]-σ_ωy, x_0_6[5]+σ_ωy)
D6 = Uniform(x_0_6[6]-σ_ωz, x_0_6[6]+σ_ωz)
x1 = zeros(1, M)
x2 = zeros(1, M)
x3 = zeros(1, M)
x4 = zeros(1, M)
x5 = zeros(1, M)
x6 = zeros(1, M)
rand!(D1, x1)
rand!(D2, x2)
rand!(D3, x3)
rand!(D4, x4)
rand!(D5, x5)
rand!(D6, x6)
x = vcat(x1, x2, x3, x4, x5, x6)
