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
    #neglect the effects of oblateness and the transport terms in ω^2
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

    if  t <= 50.0
        α = 50
    else
        α = -t+110
    end

    α = floor(Int, α)

    C_D = table_CD[α]
    C_L = table_CL[α]

    m = 600.0 # kg
    Re = 3389.5*1e3 #m
    h = r - Re #m
    D = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_D #Drag Acceleration
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

function vinhs_full_model(t, u, p)
    #no winds
    du = zeros(length(u))
    r, θ, ϕ, v, γ, ψ = u #state vector
    #theta is basically longitude
    #phi is the latitude
    σ = p[1] #bank angle (control input)

    #parameters
    req = 3396.2*1e3
    J2 = 1.96*1e-3
    ω = 7.095*10^(-5) #rad.s-1 MARS
    μ =  4.282837*1e13 #m3.s-2

    #gravity
    g_r = -μ/(r^2)*(1-1.5*J2*((req/r)^2)*(3*(sin(ϕ)^2)-1))
    g_ϕ = -3*J2*(μ/(r^2))*((req/r)^2)*sin(ϕ)*cos(ϕ)

    #Coriolis terms
    C_γ = 2*ω*cos(ψ)*cos(ϕ) #Coriolis
    C_ψ = 2*ω*(tan(γ)*sin(ψ)*cos(ϕ)-sin(ϕ)) #Coriolis

    #Transport terms
    Γ_v = r*(ω^2)*cos(ϕ)*(sin(γ)*cos(ϕ)-cos(γ)*sin(ψ)*sin(ϕ))
    Γ_γ = (r/v)*(ω^2)*cos(ϕ)*(sin(γ)*sin(ψ)*sin(ϕ)+cos(γ)*cos(ϕ))
    Γ_ψ = -(r/v)*(ω^2)*(cos(ψ)*sin(ϕ)*cos(ϕ))/(cos(γ))


    if  t <= 50.0
        α = 50
    else
        α = -t+110
    end


    α = floor(Int, α)

    C_D = table_CD[α]
    C_L = table_CL[α]

    m = 600.0 # kg
    Re = 3389.5*1e3 #m
    h = r - Re #m
    D = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_D #Drag Acceleration
    L = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_L

    #1-3 are same as simplified version
    du[1] = v*sin(γ)
    du[2] = (v/r)*(cos(γ)*cos(ψ)/cos(ϕ))
    du[3] = (v/r)*cos(γ)*sin(ψ)
    du[4] = -D+g_r*sin(γ)+g_ϕ*cos(γ)*sin(ψ) + Γ_v
    du[5] = (1/v)*(L*cos(σ)+g_r*cos(γ)-g_ϕ*sin(γ)*sin(ψ))+(v/r)*cos(γ) + C_γ + Γ_γ
    du[6] = -(1/(v*cos(γ)))*(L*sin(σ)-g_ϕ*cos(ψ)) -(v/r)*cos(γ)*cos(ψ)*tan(ϕ) +C_ψ + Γ_ψ
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

function specific_energy(V, r)
    μ = 4.282837*1e13 #m3.s-2
    E = 0.5*V^2-(μ/r) #specific energy
    return E
end


###############################################################################
#################### Regular Integration ######################################
###############################################################################

δ = 70.0*pi/180.0
r_min = 0.2
r_cone = 1.3
r_G = [0.0; 0.0; -0.189]
A_ref = pi*r_cone^2
table_CD, table_CL = drag_lift_table(δ, r_min, r_cone, r_G)

u0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
t_sim4, Z4 = rk4(vinhs_full_model, u0, [45*pi/180], 0.01, [0.0; 100.0])

Plots.plot(Z4[4, :]*1e-3, (Z4[1, :].-(3389.5*1e3))*1e-3)
Plots.plot(t_sim4, Z4[1, :])
Plots.plot(t_sim4, Z4[2, :])
Plots.plot(t_sim4, Z4[3, :])
Plots.plot(t_sim4, Z4[4, :])
Plots.plot(t_sim4, Z4[5, :])
Plots.plot(t_sim4, Z4[6, :])

#Specific and normalized energy
E = [specific_energy(Z4[4, i], Z4[1, i]) for i=1:1:length(t_sim4)]
Plots.plot(t_sim4, [specific_energy(Z4[4, i], Z4[1, i]) for i=1:1:length(t_sim4)])

Ei = specific_energy(Z4[4, 1], Z4[1, 1])
Ef = specific_energy(Z4[4, end], Z4[1, end])
Plots.plot(t_sim4, (E.-Ei)./(Ef-Ei))

a=1

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
        t_sim, Z = rk4(vinhs_model, X[:, i], [45.0*pi/180], 0.01, [t, t+dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
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
    u =[45.0*pi/180]
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

x0 = [(125+3389.5)*1e3/(3389.5*1e3); 0.0; 0.0; 7.032*1e3/(1e3*7.00); -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0/(1e3*3389.5))^2; (0.00017)^2; 0.00017^2; (1.0/(1e3*7.00))^2; (0.001)^2; (0.001)^2])
Q0 = Matrix(Q0)
A0 = inv(sqrt(Q0))
b0 = -A0*x0
t_start = 0.0
t_end = 80.0
dtt = 0.1

a=1

Alist, blist, centerlist, XX, T = uncertainty_propagation(A0, b0, t_start, t_end, dtt)

Plots.plot!(t_sim4, Z4[6, :])
Plots.scatter(T, centerlist[1, :]*1e3*3389.5)
Plots.scatter(T, centerlist[4, :]*7.00*1e3)
Plots.scatter(T, centerlist[6, :])

X_lims(XX[:, :, end])

W = inv(Alist[:, :, end])
F = eigen(W)
F.values

uncertainty(Alist)

Plots.scatter(centerlist[4, :]*7.00*1e3, centerlist[1, :]*1e3*3389.5)


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


function prop_MC_entry(X_samples, t_start, t_end, dt, p, model)
    n, M = size(X_samples)
    #saveAT = 1.0
    TT = t_start:dt:t_end
    traj = zeros(n, length(TT), M)
    for i=1:1:M
        @show(i)
        x_ini = X_samples[:, i]
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(model, x_ini, p, dt, [t_start, t_end])
        traj[:, :, i] = Z
    end
    return traj, TT
end

x0 = [(122+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0^2); (0.00017)^2; 0.00017^2; (1.0)^2; (0.001)^2; (0.001)^2])
X_samples = generate_samples(x0, Q0, 100000)
p = [45.0*pi/180]
traj, TT = prop_MC_entry(X_samples, 0.0, 90.0, 0.01, p, vinhs_model)

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

avg, var = mean_var_MC(traj)

Plots.plot!(TT, avg[6, :])
Plots.scatter(T, centerlist[6, :])
Plots.plot(TT, var[6, 6, :].^0.5)


F = eigen(inv(Alist[:, :, end]))
F.values


X_lims(traj[:, end, :])
