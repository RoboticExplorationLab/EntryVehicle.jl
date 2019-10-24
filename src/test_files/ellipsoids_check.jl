#test ellipsoid lower dimensional system
using LinearAlgebra
using DifferentialEquations
using Plots
using Convex
using SCS
using Distributions
using Random
pyplot()

function sys!(du,u,p,t)
 du[1] = u[2]
 du[2] = -u[1] #+u[1]^2
end

u0 = [0.1;0.2]
tspan = (0.0,10.0)
prob = ODEProblem(sys!,u0,tspan)
sol = solve(prob, reltol=1e-10,abstol=1e-10)

plot!(sol, vars=(1,2))


#initialization
u0 = [0.6;0.2]
Q0 = Matrix(Diagonal([0.01, 0.01]))
A0 = inv(sqrt(Q0)) #matrix at t=0.0
b0 = -A0*u0 #center at t=0.0

Δt = 3.0 #length simulation
dt = 0.05
T = 0.0:dt:Δt
n = 2

function ellipse2points(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    points = zeros(n, 6*n+1)
    M = -inv(A)*b #center of the ellipse considered
    #L = A'*A
    #F = eigen(L)
    #W = F.vectors
    #z = F.values
    W = inv(A)
    #@show(z)
    #@show(b)
    for i = 1:n
        points[:, 6*i-5] = M + 0.5*W[:, i]
        points[:, 6*i-4] = M - 0.5*W[:, i]
        points[:, 6*i-3] = M - 0.8*W[:, i]
        points[:, 6*i-2] = M + 0.8*W[:, i]
        points[:, 6*i-1] = M - W[:, i]
        points[:, 6*i] = M + W[:, i]
    end
    points[:, 6*n+1] = M
    return points #return 2n+1 points
end

function points2ellipse(X)
    #X is the set of points propagated through the non linear dynamics
    n, m = size(X);
    #Define Convex Optimization Problem using Convex.jl
    A = Semidefinite(n) #I think should be positivedefinite
    b = Variable(n)
    problem = maximize(logdet(A), [norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m])#, [A[k, j]==A[j, k] for k=1:n for j=1:n])) #[A\Matrix{Float64}(I, n, n)]))
    Convex.solve!(problem, SCSSolver(verbose = true))
    b = b.value
    A = A.value
    return A, b
end

function prop_points_continuous(X, dt)
    tspan = (0.0,dt) #autonomous system so don't care
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        prob = ODEProblem(sys!, X[:, i], tspan)
        sol = solve(prob, reltol=1e-10,abstol=1e-10)
        Xnew[:, i] = sol.u[end]
    end
    return Xnew
end


function prop(A1, b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n, 6*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1) #return a set of points
        X2 = prop_points_continuous(X1, dt) #prop_points(X1, t, dt, u, w)
        A2, b2 = points2ellipse(X2)
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


Alist, blist, centerlist, XX = prop(A0, b0)

#plot to analyse results
function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

plot_traj_center(centerlist)
j = 500
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
plot!(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!(XX[1, :, j], XX[2, :, j])
#title!("simple pendulum - no damping")
#xlabel!("x")
#ylabel!("ẋ")
savefig("pendulum - 3 sec - 6 points - 0.05 step")

#test ellipse fitting in 2D
D1 = Uniform(-3, 10)
D2 = Uniform(1, 5)
x1 = zeros(1, 5)
x2 = zeros(1, 5)
rand!(D1, x1)
rand!(D2, x2)
x = vcat(x1, x2)

n, m = size(x);

A = Semidefinite(n)
b = Variable(n)
problem = maximize(logdet(A), vcat([norm(A*x[:, i]+b, 2)<=1 for i = 1:1:m], [A[k, j]==A[j, k] for k=1:n for j=1:n]))

Convex.solve!(problem, SCSSolver(max_iters=10000))

problem.status
problem.optval

b = b.value
A = A.value
@show(A)
@show(inv(A))
@show(b)
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end

ellipse  = A \ B
scatter(x[1,:], x[2,:])
plot!(ellipse[1, :], ellipse[2, :])
Y = A\[cos(pi/2)-b[1], sin(pi/2)-b[2]]
scatter!([Y[1]], [Y[2]])
xlims!(-2.0, 11.0)
ylims!(-2.0, 11.0)

#A = [0.2135 -0.0628927; -0.0628927 0.2135]

#D = inv(A'*A)
#V = eigen(inv(A)).values


xc = -inv(A)*b
scatter!([xc[1]], [xc[2]])

scatter!([xc[1]+inv(A)[1, 1]], [xc[2]+inv(A)[2, 1]]) #corresponds to angle = 0.0 in ellipse plotting above in fact
scatter!([xc[1]-inv(A)[1, 1]], [xc[2]-inv(A)[2, 1]])
scatter!([xc[1]-inv(A)[2, 1]], [xc[2]-inv(A)[2, 2]])
scatter!([xc[1]+inv(A)[2, 1]], [xc[2]+inv(A)[2, 2]])

#Deformation of ellispes with non linear dynamics TESTS

function sys!(du,u,p,t)
 du[1] = u[2]
 du[2] = -u[1] +u[1]^2
end

#initialization
u0 = [0.1;0.2]
Q0 = Matrix(Diagonal([0.1, 0.1]))
A0 = inv(sqrt(Q0)) #matrix at t=0.0
b0 = -A0*u0 #center at t=0.0

Δt = 20.0 #length simulation
dt = 0.05
T = 0.0:dt:Δt
n = 2

function ellipse2points2(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    angles = 0.0:0.2:2*pi
    B = zeros(2, length(angles))
    for i = 1:1:length(angles)
          B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
    end
    ellipse  = A \ B
    return ellipse
end

function points2ellipse(X)
    #X is the set of points propagated through the non linear dynamics
    n, m = size(X);
    #Define Convex Optimization Problem using Convex.jl
    A = Semidefinite(n) #I think should be positivedefinite
    b = Variable(n)
    problem = maximize(logdet(A), [norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m])#, [A[k, j]==A[j, k] for k=1:n for j=1:n])) #[A\Matrix{Float64}(I, n, n)]))
    Convex.solve!(problem, SCSSolver(verbose = true))
    b = b.value
    A = A.value
    return A, b
end

function prop_points_continuous(X, dt)
    tspan = (0.0,dt) #autonomous system so don't care
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        prob = ODEProblem(sys!, X[:, i], tspan)
        sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-10,abstol=1e-10)
        Xnew[:, i] = sol.u[end]
    end
    return Xnew
end


function prop2(A1, b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    angles = 0.0:0.2:2*pi
    XX = zeros(n, length(angles), length(T))
    E = zeros(n, length(angles), length(T)) #store the propagation of the ellipses
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points2(A1, b1) #return a set of points (big set here)
        X2 = prop_points_continuous(X1, dt) #prop_points(X1, t, dt, u, w)
        A2, b2 = points2ellipse(X2)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        E[:, :, i] = X2
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        #@show(X1)
    end
    return Alist, blist, centerlist, XX, E
end

Alist, blist, centerlist, XX, E = prop2(A0, b0)

#plot to analyse results
function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

plot_traj_center(centerlist)
j = 120
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
plot!(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!(XX[1, :, j], XX[2, :, j])
scatter!(E[1, :, j], E[2, :, j])
#title!("simple pendulum - no damping")
#xlabel!("x")
#ylabel!("ẋ")
savefig("pendulum - 0.5 sec - 0.05 step")
