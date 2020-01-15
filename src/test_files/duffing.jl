#Duffing Oscillator Test case
using LinearAlgebra
using DifferentialEquations
using Plots
using Convex
using SCS
using Distributions
using Random
using Mosek, Mosek.Ext
using JuMP
using ECOS #does not solver SDPs
using MathProgBase
using Gurobi #works properly but
using MathOptInterface
using MosekTools
using LaTeXStrings
gr()

####################################
###########DUFFING OSCILLATOR#######
####################################

function duffing!(du,u,p,t)
    # p = α, β, δ, γ, ω
    α = -1.0 #p[1]
    β = 1.0 #p[2]
    δ = 0.2 #p[3]
    γ = 0.1 #p[4]
    ω = 1.0 #p[5]
    du[1] = u[2]
    du[2] =  γ*p(t)-δ*u[2]-α*u[1]-β*u[1]^3  -0.01*(u[1]-1.0)-0.01*(u[2]-0.0) #+F(t)[1]
end

function duffing(u,p,t)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = -1.0 #p[1]
    β = 1.0 #p[2]
    δ = 0.2 #p[3]
    γ = 0.1 #p[4]
    ω = 1.0 #p[5]
    du[1] = u[2]
    du[2] =  γ*cos(t)-δ*u[2]-α*u[1]-β*u[1]^3 -0.01*(u[1]-1.0)-0.01*(u[2]-0.0) #+F(t)[1]
    return du
end

function prop_points_continuous(X, dt, t)
    tspan = (t,t+dt) #autonomous system so don't care #NOT ANYMORE HERE THATS WHY
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        prob = ODEProblem(duffing!, X[:, i], tspan, M)
        sol = DifferentialEquations.solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-15,abstol=1e-15)
        Xnew[:, i] = sol.u[end]
    end
    return Xnew
end

function prop_points_rk(X, dt, t)
    tspan = [t;t+dt] #autonomous system so don't care #NOT ANYMORE HERE THATS WHY
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(duffing, X[:, i], p, 0.001, tspan) #0.00001
        Xnew[:, i] = Z[end, :]
    end
    return Xnew
end

a=1

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
        points[:, 6*i-5] = M + 0.2*W[:, i]
        points[:, 6*i-4] = M - 0.2*W[:, i]
        points[:, 6*i-3] = M - 0.4*W[:, i]
        points[:, 6*i-2] = M + 0.4*W[:, i]
        points[:, 6*i-1] = M - W[:, i]
        points[:, 6*i] = M + W[:, i]
    end
    points[:, 6*n+1] = M
    return points #return 2n+1 points
end

function ellipse2points3(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    points = zeros(n, 2*n+1)
    M = -inv(A)*b #center of the ellipse considered
    #L = A'*A
    #F = eigen(L)
    #W = F.vectors
    #z = F.values
    W = inv(A)
    #@show(z)
    #@show(b)
    for i = 1:n
        points[:, 2*i-1] = M - W[:, i]
        points[:, 2*i] = M + W[:, i]
    end
    points[:, 2*n+1] = M
    return points #return 2n+1 points
end

a = 1

function ellipse2points4(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    points = zeros(n, 2*n)
    M = -inv(A)*b #center of the ellipse considered
    #L = A'*A
    #F = eigen(L)
    #W = F.vectors
    #z = F.values
    W = inv(A)
    #@show(z)
    #@show(b)
    for i = 1:n
        points[:, 2*i-1] = M - W[:, i]
        points[:, 2*i] = M + W[:, i]
    end
    #points[:, 2*n+1] = M
    return points #return 2n+1 points
end

function ellipse2points5(A, b)
    n = length(b)
    X = zeros(n, 2*n+1)
    M = -inv(A)*b
    AA = A*A
    C =cholesky(AA)
    L = C.L
    SV = svd(L)
    U = SV.U
    Σ = Diagonal(SV.S)
    for i=1:1:n
          X[:, 2*i-1] = M-(1/Σ[i, i])*U[:, i]
          X[:, 2*i] = M+(1/Σ[i, i])*U[:, i]
    end
    X[:, 2*n+1] = M
    return X
end

function ellipse2points2(A, b)
    #A is the matrix we obtain from the previous step
    n = length(b)
    angles = 0.0:0.2:2*pi #can change that here
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
    problem = maximize(logdet(A), vcat([norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m], [A[k, j]==A[j, k] for k=1:n for j=1:n])) #[A\Matrix{Float64}(I, n, n)]))
    Convex.solve!(problem, SCSSolver(verbose = true, max_iters=50000))
    b = b.value
    A = A.value
    return A, b
end

a=1

function points2ellipse_mosek(X)
    n, m = size(X);
    s = MathOptInterface.LogDetConeTriangle(n)
    model = Model(with_optimizer(Mosek.Optimizer, MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-9), MSK_DPAR_INTPNT_CO_TOL_PFEAS=10^(-9), MSK_DPAR_INTPNT_CO_TOL_MU_RED = 10^(-10))) #MSK_DPAR_INTPNT_CO_TOL_INFEAS=10^(-12), MSK_IPAR_INTPNT_MAX_ITERATIONS=1000))
    @variable(model, A[1:n, 1:n], PSD)
    @variable(model, b[1:n])
    @variable(model , t)
    @objective(model, Max, t)
    @constraint(model, con[i = 1:m], [1.0; A*X[:, i]+b] in SecondOrderCone())
    V = [A[i, j] for j in 1:1:n for i in 1:1:j] #vectorize form of the matrix
    @constraint(model, [t;1.0;V] in s)
    @show()
    #@show(con)
    #MathOptInterface.TimeLimitSec() = 0.5
    JuMP.optimize!(model)
    a = JuMP.termination_status(model)
    #@show(getsolsta())
    @show(a)
    @show(objective_value(model))
    obj = objective_value(model)
    A = JuMP.value.(A)
    b = JuMP.value.(b)
    return A, b, obj
end

function plot_traj_center(T, centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    Plots.scatter(X, Y)
end

function uncertainty(Alist, T)
    t = length(Alist[1, 1, :])
    U = zeros(t)
    for i=1:1:t
        U[i] = tr(inv(Alist[:, :, i]))
    end
    Plots.plot(T, U, label = "Ellipses_method")
    Plots.xlabel!("time")
    Plots.ylabel!("standard_dev")
end

function rk4(f, y_0, p, dt, t_span)
    TT = t_span[1]:dt:t_span[end]
    y = zeros(length(TT), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(TT)-1
        t = TT[i]
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
    return TT, y
end

function var_points(X)
    n, m = size(X)
    X_avg = (1/m)*sum(X[:, i] for i=1:1:m)
    var = (1/m)*(sum((X[:, i]-X_avg)*(X[:, i]-X_avg)' for i=1:1:m))
    sig = sqrt(var)
    return var, sig
end

a=1

function step_elli(A1, b1, dtt, t)
    X1 = ellipse2points4(A1, b1)
    X2 = prop_points_rk(X1, dtt, t)
    A2, b2 = points2ellipse_mosek(X2)
    #@show(X2)
    #A2, b2 = DRN_algo(X2)
    return X1, X2, A2, b2#, obj
end

a = 1

function prop2(A1, b1)
    dtt = 0.01
    T = 0.0:dtt:25.0
    n = 2
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    OBJ = zeros(length(T))
    #angles = 0.0:0.2:2*pi
    #XX = zeros(n, length(angles), length(T))
    #E = zeros(n, length(angles), length(T)) #store the propagation of the ellipses
    XX = zeros(n, 4, length(T))
    E = zeros(n, 4, length(T))
    #C = zeros(length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        #=X1 = ellipse2points5(A1, b1) #return a set of points (big set here)
        #Plots.scatter!(X1[1, :], X1[2, :])
        X2 = prop_points_rk(X1, dtt, t)
        #Plots.scatter!(X2[1, :], X2[2, :])
        A2, b2, obj = points2ellipse_mosek(X2)
        angles = 0.0:0.01:2*pi
        B = zeros(2, length(angles))
        for k = 1:1:length(angles)
            B[:, k] = [cos(angles[i]) - b2[1], sin(angles[i]) - b2[2]]
        end
        ellipse  = A2 \ B
        Plots.plot!(ellipse[1, :], ellipse[2, :]) =#
        @show(i)
        X1, X2, A2, b2= step_elli(A1, b1, dtt, t)
        #if obj <= 7.0 && t >=5.0
        #    return Alist, blist, centerlist, XX, E, OBJ, T
        #end
        #@show(X1)
        #@show(X2)
        #@show(A1)
        #@show(b1)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        #OBJ[i] = obj
        E[:, :, i] = X2
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        #@show(X1)
    end
    return Alist, blist, centerlist, XX, E, OBJ, T
end

μ = [0.0] #average of the perturbation (Gaussian)
Σ = [0.01] #covariance matrix
D = MvNormal(μ, Σ)
τ = zeros(length(μ))
rand!(D, τ)
F(t) = rand!(D, τ)

#PHASE PORTRAIT AND STUFF
u0 = [0.1; 0.1]
tspan = (0.0,25.0)
p = [-1.0; 1.0; 0.2; 0.1; 1.0]
M  = t->cos(t)
prob = ODEProblem(duffing!,u0,tspan, M)
sol = DifferentialEquations.solve(prob)
Plots.plot(sol, vars=(1,2), legend = false)

Plots.plot(sol, vars=(1))
Plots.plot!(sol, vars=(2))
t_sim, Z = rk4(duffing, [0.1; 0.1], p, 0.01, [0.0;25.0])
Plots.plot(Z[:, 1], Z[:, 2])
xlabel!(L"\$\\alpha\$")
#savefig("Duffing")

#Initialization
u0 = [0.1; 0.1]
Q0 = Matrix(Diagonal([0.01^2, 0.01^2]))
A0 = inv(sqrt(Q0)) #matrix at t=0.0
b0 = -A0*u0 #center at t=0.0

#=
Δt = 50.0 #length simulation
dtt = 1.0
T = 0.0:dtt:Δt
n = 2
T = vcat(0.0:1.0:2.0, 2.1:0.1:50) =#

Alist, blist, centerlist, XX, E, OBJ, T = prop2(A0, b0)
Alist7, blist7, centerlist7, XX7, E7, OBJ7, T7 = prop2(A0, b0)

Alist8, blist8, centerlist8, XX8, E8, OBJ8, T8 = prop2(A0, b0)

Alist[:, :, 7] - Alist[:, :, 7]'

Plots.scatter(centerlist[1, :], centerlist[2, :], label = "Ellipse center", markersize = 2, legend = false)
Plots.plot!(Z[:, 1], Z[:, 2], label = "Nominal trajectory", linewidth = 1.3)

j = 2500
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    #B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
    B[:, i] = [cos(angles[i]), sin(angles[i])]
end
ellipse  = (Alist[1:2, 1:2, j])*0.2 \ B
Plots.scatter!([centerlist[1, j]], [centerlist[2, j]], markercolor = :green)
Plots.plot!(centerlist[1, j].+ ellipse[1, :],centerlist[2, j].+ ellipse[2, :], linecolor =:green, linewidth = 2.0)
Plots.scatter!([Z[floor(Int, j), 1]], [Z[floor(Int, j), 2]])
xlabel!("position")
ylabel!("velocity")
savefig("Duffing_1")


plot_traj_center(T, centerlist)
uncertainty(Alist, T)
savefig("U_0.5sec")
Plots.scatter(OBJ)
Plots.savefig("freq_comput")

a=1

anim = @animate for j=1:1:100
    angles = 0.0:0.01:2*pi
    B = zeros(2, length(angles))
    for i = 1:1:length(angles)
        B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
        #B[:, i] = [cos(angles[i]), sin(angles[i])]
    end
    ellipse  = Alist[1:2, 1:2, j] \ B
    Plots.plot(sol, vars=(1,2), legend = false)
    scatter!([centerlist[1, j]],[centerlist[2, j]] )
    Plots.plot!(ellipse[1, :], ellipse[2, :], legend = false)
    #xlims!(-0.01, 0.01)
    #ylims!(-0.01, 0.01)
    scatter!(XX[1, :, j], XX[2, :, j])
    #plot_traj_center(centerlist)
    #scatter!(E[1, :, j], E[2, :, j])
    #xlabel!("position")
    #ylabel!("velocity")
    title!("Ellipse propagation step=$(j), OBJ=$(OBJ[j])")
end
gif(anim, "test_duffing_essai.gif", fps = 1)

#Ellipse fitting at different steps
anim = @animate for j=1:1:200
    scatter(E[1, :, j], E[2, :, j], legend = false)
    angles = 0.0:0.01:2*pi
    B = zeros(2, length(angles))
    for i = 1:1:length(angles)
        B[:, i] = [cos(angles[i]) - blist[1, j+1], sin(angles[i]) - blist[2, j+1]]
    end
    ellipse  = Alist[1:2, 1:2, j+1] \ B
    plot!(ellipse[1, :], ellipse[2, :])
    scatter!(XX[1, :, j+1], XX[2, :, j+1])
    title!("Ellipse fitting t=$(j), OBJ=$(OBJ[j])")
    xlabel!("X")
    ylabel!("̇X")
end
gif(anim, "ellipse_fit.gif", fps = 1)

####################################
#######MONTE CARLO SIMULATION#######
####################################

using Distributions
using Random

u0 = [0.1; 0.1]
σ_x = 0.01 #actual value of one branch of the semi-axe
σ_y = 0.01
Q0 = Matrix(Diagonal([σ_x^2, σ_y^2]))
A0 = inv(sqrt(Q0)) #matrix at t=0.0
b0 = -A0*u0 #center at t=0.0

Δt = 25.0 #length simulation
dt = 0.01
T = 0.0:dt:Δt
n = 2

D1 = Uniform(u0[1]-σ_x/2, u0[1]+σ_x/2)
D2 = Uniform(u0[2]-σ_y/2, u0[2]+σ_y/2)
D1 = Normal(0.1, 0.003)
D2 = Normal(0.1, 0.003)
x1 = zeros(1, 1000)
x2 = zeros(1, 1000)
rand!(D1, x1)
rand!(D2, x2)
x = vcat(x1, x2)

x = generate_samples2(u0, Q0, 1000)

function prop_MC(x)
    n, m = size(x)
    #saveAT = 1.0
    tspan = (0.0,20.0)
    traj = zeros(n, 2501, m)
    for i=1:1:m
        Z = zeros(n, 2501)
        u0 = x[:, i]
        #p = [-1.0; 1.0; 0.2; 0.1; 1.0]
        t_sim, Z = rk4(duffing, u0, p, 0.01, [0.0;25.0])
        #sol= DifferentialEquations.solve(prob)#saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        traj[:, :, i] = Z'
    end
    return traj
end

function mean_var_MC(traj)
    n, t, M = size(traj)
    avg = zeros(n, t)
    var = zeros(n, n, t)
    for i=1:1:t
        @show(i)
        S = zeros(n)
        V = zeros(n, n)
        for j=1:1:M
            S += traj[:, i, j]
        end
        avg[:, i] = S/M
        for k =1:1:M
            V+= (traj[:, i, k]-avg[:, i])*(traj[:, i, k]-avg[:, i])'
        end
        var[:, :, i] = V/M
    end
    return avg, var
end

traj = prop_MC(x)

avg, var = mean_var_MC(traj)

histogram(traj[1, 2000, :], bins=100, legend = false)
savefig("MC_x_10000_20sec")
histogram(traj[2, 2000, :], bins = 100, legend = false)
savefig("MC_xvel_10000_20sec")

minimum(traj[1, 1, :])
maximum(traj[1, 1, :])

Plots.plot!(avg[1, :], avg[2, :])
Plots.scatter!(traj[1, 1550, :], traj[2, 1550, :])

Plots.plot!(traj[1, :, 1000], traj[2, :, 1000])

Plots.plot(avg[1, :]-centerlist[1, :], avg[2, :]-centerlist[2, :], legend = false)
xlabel!("position difference")
ylabel!("velocity difference")
savefig("error MC_ellipse")


Plots.plot!(T, [tr(sqrt(var[:, :, i])) for i=1:1:length(T)], label="Monte Carlo")
savefig("sigma_comparison")


sqrt(10000)

Plots.scatter!(traj[1, end, :], traj[2, end, :])
Plots.scatter!(traj[1, 1, :], traj[2, 1, :])

S = [inv(Alist[:, :, i])[1, 1] for i=1:1:length(0.0:0.01:25.0)]
Plots.plot(T, centerlist[1, :]-Z[:, 1])
Plots.plot!(T, centerlist[1, :]-Z[:, 1]+S)
Plots.plot!(T, centerlist[1, :]-Z[:, 1]-S)
Plots.plot!(T, traj[1, :, 40]-Z[:, 1])
V = [sqrt(var[1, 1, i]) for i=1:1:length(T)]
Plots.plot!(T, avg[1, :]-Z[:, 1]-V)

S = [inv(Alist[:, :, i])[2, 2] for i=1:1:length(0.0:0.01:25.0)]
Plots.plot(centerlist[2, :]-Z[:, 2])
Plots.plot!(centerlist[2, :]-Z[:, 2]+S)
Plots.plot!(centerlist[2, :]-Z[:, 2]-S)
Plots.plot!(traj[2, :, 35]-Z[:, 2])

Plots.plot(centerlist[1, :])
Plots.plot!(centerlist[1, :]+S)
Plots.plot!(centerlist[1, :]-S)
Plots.plot!(traj[1, :, 20])

#Ellipse MC comparison
anim = @animate for j=1:100:2500
    scatter(traj[1, j, :], traj[2, j, :], legend = false)
    scatter!(XX[1, :, j], XX[2, :, j])
    angles = 0.0:0.01:2*pi
    B7 = zeros(2, length(angles))
    for i = 1:1:length(angles)
        B7[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
    end
    ellipse7  = Alist[1:2, 1:2, j] \ B7
    plot!(ellipse7[1, :], ellipse7[2, :])
    #plot!(sol, vars=(1, 2))
    title!("Ellipse propagation vs MC step =$(j)")
    xlabel!("X")
    ylabel!("̇X")
end
gif(anim, "prop_vc_MC.gif", fps = 1)


#anim = @animate for j=1:1:50
j = 0
    #scatter(E7[1, :, j], E7[2, :, j], legend = false)
    #scatter!(E8[1, :, 2*j], E8[2, :, 2*j],  legend = false)
    angles = 0.0:0.01:2*pi
    B7 = zeros(2, length(angles))
    B8 = zeros(2, length(angles))
    for i = 1:1:length(angles)
        B7[:, i] = [cos(angles[i]) - blist7[1, j+1], sin(angles[i]) - blist7[2, j+1]]
        B8[:, i] = [cos(angles[i]) - blist8[1, 2*j+1], sin(angles[i]) - blist8[2, 2*j+1]]
    end
    ellipse7  = Alist7[1:2, 1:2, j+1] \ B7
    ellipse8 = Alist8[1:2, 1:2, 2*j+1] \B8
    plot!(ellipse7[1, :], ellipse7[2, :])
    plot!(ellipse8[1, :], ellipse8[2, :])
    #scatter!(XX[1, :, j+1], XX[2, :, j+1])
    title!("Ellipse fitting comparison t=$(j)")
    xlabel!("X")
    ylabel!("̇X")
#end
#gif(anim, "comparison_1_0.5.gif", fps = 1)

Alist7[:, :, 2] - Alist8[:,:, 1]


angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, 3], sin(angles[i]) - blist[2, 3]]
end
ellipse = Alist[1:2, 1:2, 3] \ B
Plots.plot!(ellipse[1, :], ellipse[2, :])
Plots.scatter!(E[1, :, ], E[2, :, end-5])
Plots.scatter!(XX[1, :, end-4], XX[2, :, end-4])


#test Mosek Break code 3943:
function step_elli(A1, b1, dtt, t)
    X1 = ellipse2points4(A1, b1)
    X2 = prop_points_rk(X1, dtt, t)
    A2, b2 = points2ellipse_mosek(X2)
    #A2, b2 = DRN_algo(X2)
    return X1, X2, A2, b2 #, obj
end

X1 = [-3329.42 -899.018 -2134.7 -2093.74; 1.1891e7 1.1891e7 1.18432e7 1.19387e7]
X2 = [-3218.55 859.37 303.248 -36.7637; 1.17365e7 -1.16494e7 -1.20092e7 -1.2067e7]
A1 = [0.000822914 -3.53009e-7; -3.53009e-7 2.09437e-5]
b1 = [5.93744, -249.042]
t = 553.09

dtt = 0.01

X1 = [1.27121 1.27389 1.28128 1.26382; 0.0221133 0.00465677 -0.0640024 0.0907725]
X2 = [1.27143 1.27394 1.28063 1.26472; 0.0218184 0.00496109 -0.0645815 0.0901944]
A1 = [2808.75 316.789; 316.789 48.6515]
b1 = [-3578.5, -403.78]
t = 43.21

dtt = 0.01

X1_af = ellipse2points4(A1, b1)
X2_af = prop_points_rk(X1_af, dtt, t)
A2, b2 = points2ellipse_mosek(X2_af)

A1 = A2
b1 = b2

################################################################################
#######################Linear Covariance Analysis###############################
################################################################################

function matrixx(x0, p)
    #x0 linearization point
    α, β, δ, γ, ω = p
    A = [0.0 1.0;
         -α-3*β*x0[1]^2-0.01 -δ-0.01]
    return A
end

function prop_1(x1, t, p, x0)
    #x0 nominal traj
    α, β, δ, γ, ω = p
    A = matrixx(x0, p)
    x2 = (A*0.01)*x1 +x1 #+[0.0; γ*cos(ω*t)]*0.01 +x1
    return x2
end

function prop(T, x0, Z)
    T = 0.0:0.01:25.0
    X = zeros(2, length(T)+1)
    x1 = x0
    X[:, 1] = x0
    for i = 1:1:length(T)
        t = T[i]
        x2 = prop_1(x1, t, p, Z[i, :])
        x1 = x2
        X[:, i+1] = x2
    end
    return X
end

x0 = [0.01; 0.01]
X = prop(T, x0, Z)

Plots.plot(X[1, :], X[2, :])
Plots.plot!(Z[:, 1], Z[:, 2])

μ0 = [0.1; 0.1]
Σ0 = [0.01^2 0.0; 0.0 0.01^2]

function lincov(μ0, Σ0)
    T = 0.0:0.01:25.0
    Q = 1e-10*[1.0 0.0; 0.0 1.0]
    ΣΣ = zeros(2, 2, length(T)+1)
    μ = zeros(2, length(T)+1)
    μ[:, 1] = μ0
    ΣΣ[:, :, 1] = Σ0
    μ1 = μ0
    Σ1 = Σ0
    for i = 1:1:length(T)
        t = T[i]
        w = rand!(MersenneTwister(1234), MultivariateNormal([0.0;0.0], Q), zeros(2))
        μ2 = (μ1 + duffing(μ1, p, t)*0.01) + w
        A = matrixx(μ1, p)
        Σ2 = (A*Σ1*A') + Q
        Σ1 = Σ2
        μ1 = μ2
        μ[:, i+1] = μ2
        ΣΣ[:, :, i+1] = Σ2
    end
    return μ, ΣΣ
end

T = 0.0:0.01:25.0
μ, ΣΣ = lincov(μ0, Σ0)
Plots.plot(μ[1, :], μ[2, :])
Plots.plot(T, [tr(ΣΣ[:, :, i])) for i=1:1:length(T)])

sqrt(ΣΣ[:, :, 200])


################################################################################
#######################Polynomial Chaos Expansion###############################
################################################################################
