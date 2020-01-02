#Vinh's model for entry vehicle
#simplified version considering:
#-The oblateness of the planet is very small
#-Transport terms due to rotation between frames are small
#STANDARD SIMPLIFIED EQUATIONS OF MOTION
#See Advances in Spacecraft Atmospheric Entry Guidance (paper)

using LinearAlgebra
using Plots
using Mosek
using MosekTools
using JuMP
using MathOptInterface

include("aero.jl")
include("quaternions.jl")

a = 1

function vinhs_model(t, u, p)
    du = zeros(length(u))
    r, θ, ϕ, v, γ, ψ = u #state vector
    #theta is basically longitude
    #phi is the latitude
    σ = p[1] #bank angle (control input)


    #@show(t)

    ω = 7.095*10^(-5) #rad.s-1 MARS
    μ =  4.282837*1e13 #m3.s-2

    g = μ/(r^2)
    C_γ = 2*ω*cos(ψ)*cos(ϕ) #Coriolis
    C_ψ = 2*ω*(tan(γ)*sin(ψ)*cos(ϕ)-sin(ϕ)) #Coriolis

    C_D = table_CD[16]
    C_L = table_CL[16]

    m = 600.0 # kg
    Re = 3389.5*1e3 #m
    h = r - Re #m
    D = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_D
    L = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_L

    du[1] = v*sin(γ)
    du[2] = (v/r)*(cos(γ)*cos(ψ)/cos(ϕ))
    du[3] = (v/r)*cos(γ)*sin(ψ)
    du[4] = -D-g*sin(γ)
    du[5] = (1/v)*(L*cos(σ)-(g-(v^2)/r)*cos(γ)) + C_γ
    du[6] = -(1/(v*cos(γ)))*(L*sin(σ)+(v^2/r)*(cos(γ)^2)*cos(ψ)*tan(ϕ))+C_ψ
    if h>0.0
        return du
    else
        return zeros(6)
    end
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
        k1 = f(t, y_star, p)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, p)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, p)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, p)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y'
end

###############################################################################
#################### Regular Integration ######################################
###############################################################################


δ = 70.0*pi/180.0
r_min = 0.2
r_cone = 0.3
r_G = [0.2; 0.0; 0.3]
A_ref = pi*r_cone^2
table_CD, table_CL = drag_lift_table(δ, r_min, r_cone, r_G)

u0 = [(122+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
t_sim4, Z4 = rk4(vinhs_model, u0, [45.0*pi/180], 0.01, [0.0; 90.0])

Plots.plot(t_sim4, Z4[1, :])
Plots.plot(t_sim4, Z4[2, :])
Plots.plot(t_sim4, Z4[3, :])
Plots.plot(t_sim4, Z4[4, :])
Plots.plot!(t_sim4, Z4[5, :])
Plots.plot!(t_sim4, Z4[6, :])


###############################################################################
######################### Uncertainty Propagation #############################
###############################################################################

function ellipse2points(A, b)
    n = length(b)
    points2 = zeros(n, 2*n+1)
    M = -inv(A)*b
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    points2[:, 2*n+1] = M
    return points2
end

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


function prop_points_rk(X, t, dt, u)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(vinhs_model, X[:, i], [45.0*pi/180], 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end

a =1

function scale_down(X2)
    n, m = size(X2)
    X3 = zeros(n, m)
    for i=1:1:m
        X3[1, i] = X2[1, i]/(1e3*3389.5)
        X3[4, i] = X2[4, i]/(1e3*7.00) #initial velocity
        X3[2:3, i] = X2[2:3, i]
        X3[5:6, i] = X2[5:6, i]
    end
    return X3
end

a=1

function scale_up(X3)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[1, i] = X3[1, i]*(1e3*3389.5)
        X4[4, i] = X3[4, i]*(1e3*7.00) #initial velocity
        X4[2:3, i] = X3[2:3, i]
        X4[5:6, i] = X3[5:6, i]
    end
    return X4
end

a=1

function uncertainty_propagation(A0, b0, t_start, t_end, dtt)
    T = t_start:dtt:t_end
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A0, b0)
        X2 = scale_up(X1)
        X3 = prop_points_rk(X2, t, dtt, u)
        X4 = scale_down(X3)
        A2, b2 = points2ellipse_mosek(X4)
        blist[:, i] = b0
        Alist[:, :, i] = A0
        centerlist[:, i] = -inv(A0)*b0
        A0 = A2
        b0 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX, T
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


function X_lims(X)
    n, m = size(X)
    lim = zeros(6, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
end

x0 = [(122+3389.5)*1e3/(3389.5*1e3); 0.0; 0.0; 7.032*1e3/(1e3*7.00); -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0/(1e3*3389.5))^2; (0.00017)^2; 0.00017^2; (1.0/(1e3*7.00))^2; (0.001)^2; (0.001)^2])
Q0 = Matrix(Q0)
A0 = inv(sqrt(Q0))
b0 = -A0*x0
t_start = 0.0
t_end = 80.0
dtt = 0.1

a=1

Alist, blist, centerlist, XX, T = uncertainty_propagation(A0, b0, t_start, t_end, dtt)

Plots.scatter(T, centerlist[1, :]*1e3*3389.5)
Plots.scatter(T, centerlist[4, :]*7.00*1e3)
Plots.scatter(T, centerlist[3, :])

X_lims(XX[:, :, end])

W = inv(Alist[:, :, end])
F = eigen(W)
F.values

uncertainty(Alist)


###############################################################################
######################### Monte Carlo Simulation ##############################
###############################################################################

using Distributions
using Random

function generate_samples(x_0, Q, M)
    #M number of samples
    n = length(x_0)
    univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    D = Product(univariate_D_vector)
    X_samples = zeros(n, M)
    rand!(D, X_samples)
    return X_samples
end

model = "vinhs_model"

function prop_MC_entry(X_samples, t_start, t_end, dt, p, model)
    n, M = size(X_samples)
    #saveAT = 1.0
    T = t_start:dt:t_end
    traj = zeros(n, length(T), M)
    for i=1:1:M
        @show(i)
        x_ini = X_samples[:, i]
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(model, x_ini, p, dt, [t_start, t_end])
        traj[:, :, i] = Z
    end
    return traj
end

X_samples = generate_samples(x0, Q0, 10)
p = [45.0*pi/180]
traj = prop_MC_entry(X_samples, 0.0, 80.0, 0.01, p)
