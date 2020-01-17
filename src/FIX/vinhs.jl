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

p = [45*pi/180]
t = 0.0
a = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]

b= 1.0

function vinhs_full_model(t, u, p)
    #no winds
    du = zeros(length(u))
    #du = zeros(eltype(u), length(u))
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
        α = 50.0
    else
        α = -t+110
    end

    #α = exp(-0.1*t)*cos(t) +15
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
    return du
end

f(u) = vinhs_full_model(0.0, u, [45.0*pi/180])

A = ForwardDiff.jacobian(f, a)


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

Plots.plot(T,  [exp(-0.1*t)*cos(t) for t in T])

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
u1 = [(80+3389.5)*1e3; 0.0; 0.0; 7800; -0.017; 1.571]
t_sim4, Z4 = rk4(vinhs_full_model, u0, [45*pi/180], 0.0001, [0.0; 80.0])

t_sim4, Z4 = rk4(vinhs_model, u1, [0.0], 0.01, [0.0;80.0])

Plots.plot(Z4[4, :]*1e-3, (Z4[1, :].-(3389.5*1e3))*1e-3)
Plots.plot(t_sim4, Z4[1, :])
Plots.plot!(t_sim4[6000:8000], Z4[2, 6000:8000])
Plots.plot(t_sim4, Z4[2, :])
Plots.plot(t_sim4, Z4[3, :])
Plots.plot(t_sim4, Z4[4, :])
Plots.plot(t_sim4, Z4[5, :])

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

function ellipse2points_s(A, b)
    n = length(b)
    points2 = zeros(n, 2*n+1)
    M = -inv(A)*b
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:n
        σ = W[:, i]/(norm[:, i])
        points2[:, 2*i-1] = M + σ
        points2[:, 2*i] = M - σ
        #@show(points2[:, 2*i])
    end
    points2[:, 2*n+1] = M
    return points2, W
end


function ellipse2points6(A, b)
    n = length(b)
    points2 = zeros(n, 6*n+1)
    M = -inv(A)*b
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:n
        points2[:, 6*i-1] = M + W[:, i]
        points2[:, 6*i-2] = M - W[:, i]
        points2[:, 6*i-3] = M + 0.5*W[:, i]
        points2[:, 6*i-4] = M - 0.5*W[:, i]
        points2[:, 6*i-5] = M + 0.8*W[:, i]
        points2[:, 6*i] = M - 0.8*W[:, i]
        #@show(points2[:, 2*i])
    end
    points2[:, 6*n+1] = M
    return points2
end

function points2ellipse_mosek(X)
    n, m = size(X);
    s = MathOptInterface.LogDetConeTriangle(n)
    model = Model(with_optimizer(Mosek.Optimizer, MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-20), MSK_DPAR_INTPNT_CO_TOL_PFEAS=10^(-20), MSK_DPAR_INTPNT_CO_TOL_MU_RED = 10^(-20)))
    #model = Model(with_optimizer(Mosek.Optimizer))
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
        t_sim, Z = rk4(vinhs_full_model, X[:, i], [45.0*pi/180], 0.01, [t, t+dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end

a =1

function scale_down(X2, S)
    n, m = size(X2)
    X3 = zeros(n, m)
    for i=1:1:m
        X3[:, i] = [X2[j,i]/S[j] for j=1:1:n]
    end
    return X3
end

a=1

function scale_up(X3, S)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[:, i] = [X3[j, i]*S[j] for j=1:1:n]
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
        X2 = scale_up(X1, S)
        X3 = prop_points_rk(X2, t, dtt, u)
        X4 = scale_down(X3, S)
        A2, b2 = points2ellipse_mosek(X4)
        blist[:, i] = b0
        Alist[:, :, i] = A0
        centerlist[:, i] = [(-inv(A0)*b0)[j]*S[j] for j=1:n]
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

function ini(x0, U0, S)
    x0_s = [x0[i]/S[i] for i=1:1:length(S)]
    U0_s = [(U0[i]/S[i])^2 for i=1:1:length(S)] #contains the sigmas squared
    Q0 = Matrix(Diagonal(U0_s))
    A0 = inv(sqrt(Q0))
    b0 = -A0*x0_s
    return A0, b0
end

function center_state_plot(T, centerlist)
    p1 = Plots.scatter(T, centerlist[1, :]*1e-3, markersize=1, ylabel = "radius [km]", xlabel = "time [s]")
    p2 = Plots.scatter(T, centerlist[2, :], markersize=1, ylabel = "longitude [rad]", xlabel = "time [s]")
    p3 = Plots.scatter(T, centerlist[3, :], markersize=1, ylabel = "latitude [rad]" , xlabel = "time [s]")
    p4 = Plots.scatter(T, centerlist[4, :], markersize=1, ylabel = "velocity [m.s-1]", xlabel = "time [s]")
    p5 = Plots.scatter(T, centerlist[5, :], markersize=1, ylabel = "flight-path-angle [rad]", xlabel = "time [s]")
    p6 = Plots.scatter(T, centerlist[6, :], markersize=1, ylabel = "heading angle [rad]", xlabel = "time [s]")
    Plots.plot(p1, p2, p3, p4, p5, p6, layout = (2, 3), legend = false)
end

S = [3389.53*1e3; 1.0; 1.0; 1e3*7.00; 1.0; 1.0]
x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
U0 = [50.0; (0.0017); (0.0017); (1.0); (0.0001); (0.0001)]
A0, b0 = ini(x0, U0, S)

t_start = 0.0
t_end = 80.0
dtt = 0.1

a=1

Alist, blist, centerlist, XX, T = uncertainty_propagation(A0, b0, t_start, t_end, dtt)

center_state_plot(T, centerlist)
uncertainty(Alist)

###############################################################################
######################### Monte Carlo Simulation ##############################
###############################################################################

using Distributions
using Random

function generate_samples(x_0, Q, M)
    #M number of samples
    n = length(x_0)
    rng = MersenneTwister(1234)
    MVN = MvNormal(x_0, Q/4)
    #univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    #D = Product(univariate_D_vector)
    X_samples = zeros(n, M)
    #rand!(D, X_samples)
    rand!(MVN, X_samples)
    return X_samples
end

function generate_samples2(x_0, Q, M)
    #M number of samples
    n = length(x_0)
    rng = MersenneTwister(1234)
    #MVN = MvNormal(x_0, Q)
    X_samples = zeros(n, M)
    #univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    #D = Product(univariate_D_vector)
    MVN = MvNormal(x_0, Q/9)
    #=c = 0
    while c < M
        @show(c)
        X = rand!(D, zeros(n))
        if (X-x_0)'*inv(Q)*(X-x_0) <= 1.0
            X_samples[:, c+1] = X
            c +=1
        end
    end =#
    X_samples = rand!(MVN, X_samples)
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

x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0^2); (0.0017)^2; (0.0017)^2; (1.0)^2; (0.0001)^2; (0.0001)^2])
X_samples = generate_samples2(x0, Q0, 10000)
p = [45.0*pi/180]
traj, TT = prop_MC_entry(X_samples, 0.0, 80.0, 0.01, p, vinhs_full_model)

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


#Plots
S = [sqrt(inv(Alist[:, :, i]*Alist[:, :, i])[1, 1]) for i=1:1:length(T)]
Plots.plot(T, (centerlist[1, :]*3389.5).-3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(T, (centerlist[1, :]-S.-1.0)*3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(T, (centerlist[1, :]+S.-1.0)*3389.5-((Z4[1, 1:100:end]*1e-3).-3389.5))
Plots.plot!(TT, ((traj[1, :, 60])-avg[1, :]), linewidth=0.1, legend = false)
V = [sqrt(var[1, 1, i]) for i=1:1:length(TT)]
Plots.plot!(TT, (avg[1, :]-Z4[1, :]-3*V))
Plots.plot!(TT, (avg[1, :]-Z4[1, :]+3*V))

T = 0.0:0.1:80.0
pa = 10
S = [inv(Alist[:, :, i])[1, 1] for i=1:1:length(T)]*3e3*1e3
S = [sqrt(inv(Alist[:, :, i]*Alist[:, :, i])[1, 1]) for i=1:1:length(T)]*3389.5*1e3
Plots.plot(T, centerlist[1, :]-Z4[1, 1:pa:end])
Plots.plot!(T, centerlist[1, :]-Z4[1, 1:pa:end]-S)
Plots.plot!(T, centerlist[1, :]-Z4[1, 1:pa:end]+S)
V = [sqrt(var[1, 1, i]) for i=1:1:length(TT)]
Plots.plot!(TT, avg[1, :]-Z4[1, :]-3*V, linestyle = :dot, color = :black)
Plots.plot!(TT, avg[1, :]-Z4[1, :]+3*V, linestyle = :dot, color = :black, legend = false)
Plots.plot(TT, [traj[1, :, i]-Z4[1, :] for i=1:3000], linewidth = 0.01, color = :blue)

Plots.plot!(T,-S, linestyle = :dash, color = :red)
Plots.plot!(T,+S, linestyle = :dash, color = :red)
Plots.plot!(TT, -3*V, linestyle = :dot, color = :black)
Plots.plot!(TT, +3*V, linestyle = :dot, color = :black, legend = false)


S = [sqrt(inv(Alist[:, :, i]*Alist[:, :, i])[4, 4]) for i=1:1:length(T)]*1e3*7.00
Plots.plot(T, centerlist[4, :]-Z4[4, 1:pa:end])
Plots.plot(T, centerlist[4, :]-S-Z4[4,1:pa:end], linestyle = :dash, color = :red)
Plots.plot!(T, (centerlist[4, :]+S-Z4[4, 1:pa:end]), linestyle = :dash, color = :red)
V = [sqrt(var[4, 4, i]) for i=1:1:length(TT)]
Plots.plot!(TT, avg[4, :]-Z4[4, :]-3*V, linestyle = :dot, color = :black)
Plots.plot!(TT, avg[4, :]-Z4[4, :]+3*V, linestyle = :dot, color = :black, legend = false)
Plots.plot!(TT, [traj[4, :, i]-Z4[4, :] for i=1:1000], linewidth = 0.1, color = :blue)

Plots.plot(T,-S, linestyle = :dash, color = :red)
Plots.plot!(T,+S, linestyle = :dash, color = :red)
Plots.plot!(TT, -3*V, linestyle = :dot, color = :black)
Plots.plot!(TT, +3*V, linestyle = :dot, color = :black)

c=3
S = [sqrt(inv(Alist[:, :, i]*Alist[:, :, i])[c, c]) for i=1:1:length(T)]*1.0
Plots.plot(T, centerlist[c, :]-Z4[c, 1:pa:end])
Plots.plot(T, centerlist[c, :]-S-Z4[c,1:pa:end], linestyle = :dash, color = :red)
Plots.plot!(T, (centerlist[c, :]+S-Z4[c, 1:pa:end]), linestyle = :dash, color = :red)
V = [sqrt(var[c, c, i]) for i=1:1:length(TT)]
Plots.plot!(TT, avg[c, :]-Z4[c, :]-3*V, linestyle = :dot, color = :black)
Plots.plot!(TT, avg[c, :]-Z4[c, :]+3*V, linestyle = :dot, color = :black, legend = false)
Plots.plot!(TT, [traj[c, :, i]-Z4[c, :] for i=1:1000], linewidth = 0.1, color = :blue)


function deviation(Alist)
    n, n, m = size(Alist)
    σ = zeros(n, m)
    for i=1:1:m
        W = inv(Alist[:, :, i])
        F = eigen(W)
        X = F.values
        σ[:, i] = X
    end
    return σ
end
σ = deviation(Alist)

###############################################################################
################################Test MC########################################

function uncertainty_propagation_MC(A0, b0, t_start, t_end, dtt, M)
    T = t_start:dtt:t_end
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n, M+13, length(T))
    u =[45.0*pi/180]
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        Q = (inv(A0))^2
        X12 = generate_samples2(-inv(A0)*b0, Q, M)
        X13 = ellipse2points(A0, b0)
        X1 = hcat(X12, X13)
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


x0 = [(125+3389.5)*1e3/(1e3*3389.5); 0.0; 0.0; 7.032*1e3/(1e3); -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0/(1e3*3389.5))^2; (0.0017)^2; (0.0017)^2; (1.0/(1e3))^2; (0.0001)^2; (0.0001)^2])
#Q0 = Diagonal([((10.0/(1e3*3389.5))^2); (0.00017)^2; (0.00017)^2; (1.0/1e3)^2; (0.017)^2; (0.017)^2])
Q0 = Matrix(Q0)
A0 = inv(sqrt(Q0))
b0 = -A0*x0
t_start = 0.0
t_end = 50.0
dtt = 0.1
M = 100

a=1

Alist, blist, centerlist, XX, T = uncertainty_propagation_MC(A0, b0, t_start, t_end, dtt, M)


################################################################################
##################################PCE###########################################
################################################################################


n = 40 #Number of rec points needed
d = 5 #higher degree of multivariate polynomials
op = OrthoPoly("gaussian",  d , Nrec=n) #if needed Nrec enables to compute more recurrence coefficients
#op = GaussOrthoPoly(d)
opq = GaussOrthoPoly(d; Nrec=n) #probabilist Hermite
N = 6 #number of random inputs
mop = MultiOrthoPoly([opq for i=1:N], d)
P = mop.dim #total number of Polynomials
mop.ind
showbasis(opq; sym="ξ")



####################################
##############Sampling##############
####################################

Ms = 2000 #number of sampled trajectories

Q0 = Matrix(Diagonal([U0[i]^2 for i=1:1:length(U0)]))
D = MvNormal(x0, Q0/9)
ξ = zeros(length(x0), Ms)
rand!(D, ξ) #one column is one state sample in x

#integration parameters
t0 = 0.0
dt = 0.1 #stored every unit (no implication on solver)
tf = 80.0
t_sim = t0:dt:tf
w = [0.0158*10^9; 0.0; 0.0; 0.0] #can be varied later

δ = 70*pi/180
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_G) #offline coefficients computation



function generate_samples_PCE(t0, tf, dt, ξ, Ms)
    samples = zeros(6, length(t_sim), Ms)
    for i=1:1:Ms
        #one column of Z contains the state at a specific time
        t_sim, Z = rk4(vinhs_full_model, ξ[:, i], [45.0*pi/180], dt, [t0; tf])
        samples[:, :, i] = Z
        @show(i)
    end
    return samples
end

samples = generate_samples_PCE(t0, tf, dt, ξ, Ms)

####################################
######Coefficients Computation######
####################################

function compute_Phi_matrix(Ms, P, mop, samples, x0, Q0)
    Phi = zeros(Ms, P)
    for i=1:Ms
        for j=1:1:P
            vec = [(ξ[k, i]-x0[k])/sqrt(Q0[k, k]) for k=1:1:6]
            res = PolyChaos.evaluate(mop.ind[j, :], vec, mop)[1]
            Phi[i, j] = res
        end
    end
    return Phi
end

function compute_coeff_PCE_step(samples, t, A) #t is the time step we want to look at t belongs to
    #function computes PCE coefficients for time step t
    C = zeros(P, 6) #contains coefficients for one time step
    for i=1:1:6
        C[:, i] = (A*samples[i, t, :])'
    end
    return C
end

function compute_coeff_PCE_full(samples, Phi)
    n, t, Ms = size(samples)
    C = zeros(P, n, t)
    A = pinv(Phi)
    T = t0:dt:tf
    for j = 1:1:length(T)
        c = compute_coeff_PCE_step(samples, j, A)
        C[:, :, j] = c
    end
    return C
end

Phi = compute_Phi_matrix(Ms, P, mop, samples, x0, Q0)
C = compute_coeff_PCE_step(samples, 80, Phi)
CC = compute_coeff_PCE_full(samples, Phi)

function mean_var_PCE(CC, T)
    t = length(T)
    m = zeros(6, t)
    var = zeros(6, t)
    for j=1:1:t
        m[:, j] = CC[1, :, j]
        for i = 1:1:6
            var[i, j] = sum(CC[k, i, j]^2 for k=2:1:P)
        end
    end
    return m, var
end

T = t0:dt:tf
m_PCE, var_PCE = mean_var_PCE(CC, T)
Plots.plot(T, m_PCE[1, :])
Plots.plot!(TT, Z4[1, :])
VV = [sqrt(var_PCE[1,i]) for i=1:1:length(T)]
Plots.plot!(T, -VV, linestyle = :dashdotdot, color = :green)
Plots.plot!(T, +VV, linestyle = :dashdotdot, color = :green)


Plots.plot(T, CC[1, 1, :])
Plots.plot(T, CC[1, 5, :])

size(samples)
size(Phi)

pinv(Phi)*samples[1, 50, :]

################################################################################
#####################Linear Covariance Analysis#################################
################################################################################

function vinhs_full_model_lincov(t, u, p)
    #no winds
    #du = zeros(length(u))
    du = zeros(eltype(u), length(u))
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
        α = 50.0
    else
        α = -t+110
    end

    #α = exp(-0.1*t)*cos(t) +15
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
    return du
end


function prop_dyn(u0, tstart, tend, dt)
    T = tstart:dt:tend
    p = [45*pi/180]
    u1 = u0
    U = zeros(length(u0), length(T)+1)
    U[:, 1] = u1
    for i=1:1:length(T)
        t = T[i]
        u2 = u1 + vinhs_full_model(t, u1, p)*dt
        U[:, i] = u2
        u1 = u2
    end
    return U
end

u0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
U = prop_lincov(u0, 0.0, 1.0, 1e-6)

Plots.plot(U[1, :])

#initial dispersion from nominal
u0
Q0 = Matrix(Diagonal(U0.^2)/9)

function lincov_prop(x0, Q0, Z4, T)
    E = zeros(6, length(T))
    Q = zeros(6, 6, length(T))
    for i=1:1:length(T)
        t = T[i]
        p = [45*pi/180]
        f(u) = vinhs_full_model_lincov(t, u, p)
        A = ForwardDiff.jacobian(f, Z4[:, i])
        Φ = exp(A*t)
        x = Φ*x0
        E[:, i] = x
        Q[:, :, i] = Φ*Q0*Φ'
    end
    return E, Q
end

T = 0.0:0.001:80.0
E, Q = lincov_prop(u0, Q0, Z4, T)

t_sim4, Z4 = rk4(vinhs_full_model, u0, [45*pi/180], 0.01, [0.0; 80.0])
Plots.plot!(T, [sqrt(Q[4, 4, i]) for i=1:1:length(T)])
Plots.plot!(T, -[sqrt(Q[4, 4, i]) for i=1:1:length(T)])


Plots.plot(T, E[1, :])
