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
        t_sim, Z = rk4(dyna_coeffon_COM_on_axis, X[:, i], u, 0.01, [0.0, dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
        #rk4(dyna_coeffoff, X[:, i], u, 0.001, [0.0, dt])
        @show(i)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end


function ellipse2points(A, b, x0)
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
    return points3 #return 2n+1 points
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


function propagation(A1, b1)
    Δt = 180.0
    dtt = 1.0
    T = 0.0:dtt:Δt
    Re = 3389.5*1e3
    x0 = [(3389.5+125)*1e3/(1e3*3389.5); 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci/(1e3*7.00); 0.0; 0.0; 0.0]
    n = length(b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1, x0) #return a set of points in dim 13
        X2 = scale_up(X1)
        X3 = prop_points_rk(X2, t, dtt, u)
        x1 = X3[:, end]
        X_12 = points13_12(X3, x1)
        X4 = scale_down(X_12)
        A2, b2 = points2ellipse_mosek(X4)
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

function scale_up(X3)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[1:3, i] = X3[1:3, i]*(1e3*3389.5)
        X4[8:10, i] = X3[8:10, i]*(1e3*7.00) #initial velocity
        X4[11:13, i] = X3[11:13, i]
        X4[4:7, i] = X3[4:7, i]
    end
    return X4
end

function scale_down(X3)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[1:3, i] = X3[1:3, i]/(1e3*3389.5)
        X4[7:9, i] = X3[7:9, i]/(1e3*7.00) #initial velocity
        X4[10:12, i] = X3[10:12, i]
        X4[4:6, i] = X3[4:6, i]
    end
    return X4
end
