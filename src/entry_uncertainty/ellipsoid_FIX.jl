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
pyplot()
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
Δt = 60 #length simulation

 ####################################
 #####Dynamics - Integration#########
 ####################################

t_sim_nom, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)

 #new part for offline aerodynamic coefficients computation
δ = 70*pi/180
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_G)

t_sim_nom, Z = integration2(dyna_coeffoff_inplace!, x0, Δt)

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
Q0 = Diagonal([(1.0/Re)^2;(1.0/Re)^2; (1.0/Re)^2; 0.01^2; 0.01^2; 0.01^2; 0.05^2; 0.05^2;0.05^2; 0.0000001; 0.0000001; 0.0000001]) #here we assume the change in a(vector part of the change in quaternion)
#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_12

Δt = 80.0 #length simulation
dt = 0.5
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
        t_sim, Z = rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
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
Q0 = Diagonal([(1.0/Re)^2;(1.0/Re)^2; (1.0/Re)^2; 0.01^2; 0.01^2; 0.01^2; 0.05^2; 0.05^2;0.05^2; 0.0000001; 0.0000001; 0.0000001]) #here we assume the change in a(vector part of the change in quaternion)

#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0_12

Δt = 100.0 #length simulation
dt = .5
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
        #@show(X1)
        X2 = prop_points_last(X1, dt, u, w)
        x1 = X2[:, end] #we store the last (previous center propagated)
        #@show(X2)
        X3 = points13_12(X2, x1)
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

Alist, blist, centerlist, XX, WW = propagation(A1, b1)

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
        U[i] = tr(Alist[:, :, i])
    end
    Plots.plot(U)
    Plots.xlabel!("time")
    Plots.ylabel!("uncertainty matrix trace")
end


plot_traj_center(centerlist)
uncertainty(Alist)

Plots.plot(centerlist[1, :])

maximum(centerlist[1, :])
minimum(centerlist[1, :])

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
