#test satellite orbit
using Plots
using LinearAlgebra
using Mosek
using JuMP
using MathOptInterface
using MosekTools

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

function atm_dens(r, ρ0)
    Re = 6378.137 #km
    H = 10.0 #km characteristic atm length for earth
    ρ = ρ0*exp(-(r-Re)/H)
    return ρ
end

function dyna_sat(t, x, u)
    #no control u in our case here
    x_dot = zeros(length(x))
    m = 1500.0 #kg
    CD = 2.3 #NA
    A = 20e-6 #km2
    μ = 398600.4418 #km3.s-2
    ρ0 = 1.225e9 #kg.km-3
    x_dot[1:3] = x[4:6]
    r = norm(x[1:3])
    v = norm(x[4:6])
    B = CD*A/m #ballistic coeff
    a_drag = -0.5*v*B*atm_dens(r, ρ0)*x[4:6]
    x_dot[4:6] = -mu*x[1:3]*(1/(r^3)) #+ a_drag
    return x_dot
end

#ini values
mu = 398600.4418
Re = 6378.137
rp = 200+Re
e = 0.0
a = rp/(1-e)
vp = (2*((-mu/(2*a))+(mu/rp)))^(0.5) #km.s-1
r_ECI = [rp; 0.0; 0.0] #[km] Starts at perigee
v_ECI = [0.0; vp; 0.0] #[km.s-1]
X0 = [r_ECI; v_ECI]

t_span = [0.0, 2000.0]
dt = 0.01
p = [0.0]
y_0 = X0
t_sim, Z = rk4(dyna_sat, y_0, p, dt, t_span)

Plots.plot!(Z[1, :], Z[2, :])
Plots.plot(t_sim, [norm(Z[1:3, i])-Re for i=1:1:length(t_sim)])


################################################################################
#########################Ellipsoid fit##########################################
################################################################################

#uncertainty init
x0 = y_0
Q0 = Diagonal([0.01^2; 0.01^2; 0.01^2; 0.0001^2; 0.0001^2; 0.0001^2])
Q0 = Matrix(Q0)
A0 = inv(sqrt(Q0))
b0 = -A0*x0

function ellipse2points(A, b)
    n = length(b)
    points2 = zeros(n, 2*n)
    M = -inv(A)*b
    W = inv(A)
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
    end
    return points2
end

ellipse2points(A0, b0)

a =1.0

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

a = 1.0

function prop_points_rk(X, dt)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    u = [0.0]
    for i=1:1:m
        t_sim, Z = rk4(dyna_sat, X[:, i], u, 0.01, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end

function propagation(A1, b1)
    dt = 2.0
    T = 0.0:dt:2000.0
    n = length(b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(6, 2*n, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1)
        X2 = prop_points_rk(X1, dt)
        A2, b2 = points2ellipse_mosek(X2)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX, T
end


Alist, blist, centerlist, XX, T = propagation(A0, b0)


#############Results#################

function plot_traj_center(T, centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]
        Y[j] = centerlist[2, j]
    end
    Plots.scatter(X, Y)
end

function X_lims(X)
    n, m = size(X)
    lim = zeros(6, 2)
    for i =1:1:n
        lim[i, 1], lim[i, 2] = minimum(X[i, :]), maximum(X[i, :])
    end
    return lim
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

plot_traj_center(T, centerlist)
Plots.scatter(T, [norm(centerlist[1:3, i])-Re for i=1:1:length(T)])
X_lims(XX[:, :, end])
uncertainty(Alist)
