# Modular version of the ellipsoid fiiting

function ellipse2points(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    n = length(b)
    points2 = zeros(n, 2*n+1)
    points3 = zeros(13, 2*n+1)
    M = -inv(A)*b
    C = zeros(13)
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    C[8:10] = M[7:9]
    C[11:13] = x0[11:13]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    for i =1:n
        points2[:, 2*i-1] = M + W[:, i]
        points2[:, 2*i] = M - W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:2*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:10, i] = points2[7:9, i]
        points3[11:13, i] = x0[11:13]
    end
    points3[:, 2*n+1] = C
    return points3
end

function points13_9(X_13, x1)
    #takes points in 13 and return points in 9
    n, m = size(X_13)
    X_9 = zeros(9, m)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    for i=1:1:m
        δq = qmult(qconj(x1[4:7]), X_13[4:7, i])
        X_9[1:3, i] = X_13[1:3, i]
        X_9[4:6, i] = 2*V*log_quat(δq)
        X_9[7:9, i] = X_13[8:10, i]
    end
    return X_9, X_13[:, end]
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


function prop_points_rk(X, t, dt, u)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(model, X[:, i], u, 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
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

#n = dimension of the ellipsoid space
#m = number of extracted and propagated points at each step
#N = dimension of the state vector associated with the dynamics
#fit = 'mosek', 'DRN', 'SCS'
#model = "duffing', 'entry_6DOF', 'entry_vinhs'
#dtt = frequency of computation of the ellipsoids
#x0 nominal initial condition (careful if reference, attitude needs to be 0.0)
#integrator (later between RK4 and DiffEq)

function uncertainty_propagation(A0, b0, t_start, t_end, dtt, x0, m)
    T = t_start:dtt:t_end
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = prop_points_rk(model, X1, t, dtt, u, w)
        x1 = X2[:, end]
        X_9, x0 = points13_9(X2, x1)
        A2, b2 = points2ellipse_mosek(X_9)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
        x0 = x1
    end
    return Alist, blist, centerlist, XX, T
end
