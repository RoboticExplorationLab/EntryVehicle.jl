#Carlemann linearization on Duffing Oscillator
using Plots

function S(i, p)
    S = 0
    for k =0:1:i-1
        S += (p+1)-k
    end
    return S
end

function map_index(i, j, p)
    if i == 0
        return (j+1)
    else
        return S(i, p) + (j+1)
    end
end

a=1

function carlemann_matrix(p, t)
    #p is the order at which we truncate the linearization (polynomial order)
    N = Int64((p+1)*(p+2)/2)
    A = zeros(N, N)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p) #number of the line we are working at
            if k <= N && k >=1
                A[k, k] = -δ*j
                k1 = map_index(i+1, j-1, p)
                if k1 <= N && k1 >=1
                    A[k, k1] = -j*α
                end
                k2 = map_index(i+3, j-1, p)
                if k2 <= N && k2 >=1
                    A[k, k2] = -j*β
                end
                k3 = map_index(i, j-1, p)
                if k3 <= N && k3 >=1
                    A[k, k3] = -j*γ*cos(t)
                end
                k4 = map_index(i-1, j+1, p)
                if k4 <= N && k4 >=1
                    A[k, k4] = i
                end
            end
        end
    end
    return A
end

α = -1.0
β = 1.0
γ = 0.1
δ = 0.2
A = carlemann_matrix(2, pi)

function basis(x, p)
    N = Int64((p+1)*(p+2)/2)
    base = zeros(N)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p)
            base[k] = (x[1]^i)*(x[2]^j)
        end
    end
    return base
end

function integration_carlemann(x_ini, t_ini, t_end, dt)
    T = t_ini:dt:t_end
    X = zeros(length(x_ini), length(T))
    X[:, 1] = x_ini
    for q=1:1:length(T)-1
        A = carlemann_matrix(p, T[q])
        Phi = exp(A*dt)
        k1 = map_index(1, 0, p)
        k2 = map_index(0, 1, p)
        X[:, q+1] = vcat(Phi[k1, :]', Phi[k2, :]')*basis(X[:, q], p)
    end
    return X
end

t_ini = 0.0
t_end = 100.0
dt = 0.01
p = 4 #order 3 is enough when dt small enough
x_ini = [0.1 10.0]
N = Int64((p+1)*(p+2)/2)
X = integration_carlemann(x_ini, t_ini, t_end, dt)
Plots.plot(X[1, :], X[2, :], legend = false)
savefig("duff_int")


σ_x1 = 0.01
σ_x2 = 0.01

n = 15

T = 0.0:1.0:10.0
Q_14 = Matrix(Diagonal(0.0001*ones(N-1)))
x_ini14 = basis(x_ini, p)[2:end]

A0 = inv(sqrt(Q_14)) #matrix at t=0.0
b0 = -A0*x_ini14 #center at t=0.0

a = 1

####################################
######## Apply Ellipsoid Fit #######
####################################


function ellipse2points4(A, b)
    #A is the matrix we obtain from the previous step
    nn = length(b)
    points = zeros(nn, 2*nn)
    M = -inv(A)*b #center of the ellipse considered
    #L = A'*A
    #F = eigen(L)
    #W = F.vectors
    #z = F.values
    W = inv(A)
    #@show(z)
    #@show(b)
    for i = 1:nn
        points[:, 2*i-1] = M - W[:, i]
        points[:, 2*i] = M + W[:, i]
    end
    #points[:, 2*n+1] = M
    return points #return 2n+1 points
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
    #@show(con)
    #MathOptInterface.TimeLimitSec() = 0.5
    JuMP.optimize!(model)
    a = JuMP.termination_status(model)
    @show(a)
    @show(objective_value(model))
    obj = objective_value(model)
    A = JuMP.value.(A)
    b = JuMP.value.(b)
    return A, b, obj
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

a=1

function prop_points_carl(X, dt, t)
    T = t:0.01:t+dt
    m = length(X[1, :])
    Xnew = zeros(size(X))
    XX = zeros(n, m, length(T))
    for q=1:1:length(T)-1
        A = carlemann_matrix(p, T[q])
        Phi = exp(A*dt)
        for i = 1:1:m
            XX[:, i, 1] = X[:, i]
            XX[:, i, q+1] = Phi*XX[:, i, q]
        end
        #X[:, q+1] = hcat(Phi[k1, :], Phi[k2, :])'*basis(X[:, q], p)
    end
    Xnew = XX[:, :, end]
    return Xnew
end

a=1

function points15_14(X)
    return X[2:end, :]
end

function points14_15(X)
    n, m = size(X)
    return vcat(ones(m)', X)
end

a = 1

function prop2(A1, b1)
    blist = zeros(n-1, length(T))
    Alist = zeros(n-1, n-1, length(T))
    centerlist = zeros(n-1, length(T))
    OBJ = zeros(length(T))
    angles = 0.0:0.2:2*pi
    #XX = zeros(n, length(angles), length(T))
    #E = zeros(n, length(angles), length(T)) #store the propagation of the ellipses
    XXX = zeros(n-1, 28, length(T))
    E = zeros(n, 28, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points4(A1, b1) #return a set of points (big set here)
        @show(n)
        @show(size(b1))
        X2 = points14_15(X1)
        @show(size(X2))
        X3 = prop_points_carl(X2, dt, t)
        X4 = points15_14(X3)
        @show(size(X4))
        A2, b2, obj = points2ellipse_mosek(X4)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        OBJ[i] = obj
        E[:, :, i] = X3
        A1 = A2
        b1 = b2
        XXX[:, :, i] = X1
        #@show(X1)
    end
    return Alist, blist, centerlist, XXX, E, OBJ
end

Alist, blist, centerlist, XX, E, OBJ = prop2(A0, b0)

########################################################
######Ellipse Propagation using Linearize dynamics######
########################################################

E_0 = Matrix(Diagonal(0.0001*ones(N))) #uncertainty in 15 dimension

function prop_uncertainty_lin(E_0)
    T = 0.0:1.0:50.0
    EE = zeros(N, N, length(T))
    EE[:, :, 1] = E_0
    for i = 1:1:length(T)-1
        t = T[i]
        A = carlemann_matrix(p, t)
        EE[:, :, i+1] = A*EE[:, :, i]*A'
    end
    return EE
end

prop_uncertainty_lin(E_0) #big numbers here, not good, Ellipsoids grow unbounded as if there was chaos

#Implement SOS to see if we can find a matrix P satisfying Lyapunov equation
#with the state vector being the one in 15 dimensions basically.

Plots.savefig("Duffing_Carlemann")
