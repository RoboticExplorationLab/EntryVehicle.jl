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
        t_sim, Z = rk4(dyna_coeffoff_COM_on_axis, X[:, i], u, 0.01, [t, t+dt])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
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
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    for i=1:m
        X_12[1:3, i] = X[1:3, i]
        δq = qmult(qconj(x1[4:7]), X[4:7, i])
        #X_12[4:6, i] = X[5:7, i] #nah need deltaq with respct to a ref quat
        #a = qmult(qconj(x1[4:7]), X[4:7, i])
        X_12[4:6] = 2*V*log_quat(δq)
        X_12[7:12, i] = X[8:13, i]
    end
    return X_12
end

function ellipse_propagation(A0, b0, tstart, tend, dtt)
    T = tstart:dtt:tend
    x0_i = [(3389.5+125)*1e3/(1e3*3389.5); 0.0; 50.0/(1e3*3389.5); Q[1]; Q[2]; Q[3]; Q[4]; v_eci/(1e3*7.00); 0.0; 0.0; 0.0]
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    ref = zeros(13, length(T))
    u = [0.0]
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A0, b0, x0_i) #return a set of points in dim 13
        X2 = scale_up(X1)
        X3 = prop_points_rk(X2, t, dtt, u)
        x1 = X3[:, end]
        X_12 = points13_12(X3, x1)
        X4 = scale_down(X_12)
        A2, b2 = points2ellipse_mosek(X4)
        blist[:, i] = b0
        Alist[:, :, i] = A0
        centerlist[:, i] = -inv(A0)*b0
        A0 = A2
        b0 = b2
        XX[:, :, i] = X1
        ref[:, i] = x0_i
        x0_i = scale_down_13(x1)
    end
    return Alist, blist, centerlist, XX, ref, T
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

function scale_down_13(X)
    #X is a single vector
    n = length(X)
    X2 = zeros(n)
    X2[1:3] = X[1:3]/(1e3*3389.5)
    X2[8:10] = X[8:10]/(1e3*7.00) #initial velocity
    X2[11:13] = X[11:13]
    X2[4:7] = X[4:7]
    return X2
end


function point12_13(X, xc)
    #X is a point in 12 dimension and scaled for fitting (ellipsoid space)
    #xc point of reference
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    #X_12_up = scale_up(X)
    X13 = zeros(13)
    X_12_up = X
    X13[1:3] = X_12_up[1:3]
    X13[4:7] = qmult(xc[4:7], exp_quat(V'*X_12_up[4:6]/2))
    X13[8:13] = X_12_up[7:12]
    return X13
end


function ellipse2points2(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points2 = zeros(n, 6*n)
    points3 = zeros(n+1, 6*n+1)
    M = -inv(A)*b
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    #C[1:4] = [1.0; 0.0; 0.0; 0.0]
    C[8:13] = M[7:12]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)
    for i =1:n
        points2[:, 6*i-5] = M + 0.5*W[:, i]
        points2[:, 6*i-4] = M - 0.5*W[:, i]
        points2[:, 6*i-3] = M - 0.8*W[:, i]
        points2[:, 6*i-2] = M + 0.8*W[:, i]
        points2[:, 6*i-1] = M - W[:, i]
        points2[:, 6*i] = M + W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:6*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:13, i] = points2[7:12, i]
    end
    points3[:, 6*n+1] = C
    return points3 #return 2n+1 points
end

function ellipse2points3(A, b, x0)
    #A is the matrix we obtain from the previous step
    #x is the center of the ellipsoid in 7 dimensions
    #y is the center of the ellipsoid in 6 dimensions
    n = length(b)
    V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    points2 = zeros(n, 10*n)
    points3 = zeros(n+1, 10*n+1)
    M = -inv(A)*b
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = qmult(x0[4:7], exp_quat(V'*M[4:6]/2)) #still need the reference
    #C[1:4] = [1.0; 0.0; 0.0; 0.0]
    C[8:13] = M[7:12]
    W = inv(A) #W is D^(0.5) if A coming from convex problem is symmetric...
    #@show(W)
    for i =1:n
        points2[:, 10*i-9] = M + 0.6*W[:, i]
        points2[:, 10*i-8] = M - 0.6*W[:, i]
        points2[:, 10*i-7] = M + 0.2*W[:, i]
        points2[:, 10*i-6] = M - 0.2*W[:, i]
        points2[:, 10*i-5] = M + 0.5*W[:, i]
        points2[:, 10*i-4] = M - 0.5*W[:, i]
        points2[:, 10*i-3] = M - 0.8*W[:, i]
        points2[:, 10*i-2] = M + 0.8*W[:, i]
        points2[:, 10*i-1] = M - W[:, i]
        points2[:, 10*i] = M + W[:, i]
        #@show(points2[:, 2*i])
    end
    for i =1:10*n
        points3[1:3, i] = points2[1:3, i]
        points3[4:7, i] = qmult(x0[4:7], exp_quat(V'*points2[4:6, i]/2))
        points3[8:13, i] = points2[7:12, i]
    end
    points3[:, 10*n+1] = C
    return points3 #return 2n+1 points
end

function uncert_mat(Alist)
    n,m,p = size(Alist)
    QQ = zeros(n, m, p)
    for j=1:1:p
        Q = inv(Alist[:, :, j]*Alist[:, :, j]')
        QQ[:, :, j] = Q
    end
    return QQ
end

################################################################################
#########################NEW METHOD - NO EULER DIFFERENCE ######################
################################################################################

function ellipse2points_new(A, b)
    n = length(b)
    points2 = zeros(n, 2*n)
    points3 = zeros(n+1, 2*n+1)
    M = -inv(A)*b
    C = zeros(13, 1) #center of the computed ellipoid in the real 13 dimension
    C[1:3] = M[1:3]
    C[4:7] = euler2quat(M[4:6])
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
        points3[4:7, i] = euler2quat(points2[4:6, i])
        points3[8:13, i] = points2[7:12, i]
    end
    points3[:, 2*n+1] = C
    return points3 #return 2n+1 points
end


function points13_12_new(X)
    n, m = size(X);
    X_12 = zeros(n-1, m)

    for i=1:m
        X_12[1:3, i] = X[1:3, i]
        X_12[4:6, i] = quat2euler(X[4:7, i])
        X_12[7:12, i] = X[8:13, i]
    end
    return X_12
end


function ellipse_propagation_new(A0, b0, tstart, tend, dtt)
    T = tstart:dtt:tend
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    ref = zeros(13, length(T))
    u = [0.0]
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points_new(A0, b0) #return a set of points in dim 13
        X2 = scale_up(X1)
        X3 = prop_points_rk(X2, t, dtt, u)
        x1 = X3[:, end]
        X_12 = points13_12_new(X3)
        X4 = scale_down(X_12)
        A2, b2 = DRN_algo(X4)
        blist[:, i] = b0
        Alist[:, :, i] = A0
        centerlist[:, i] = -inv(A0)*b0
        A0 = A2
        b0 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX, T
end

#test

v_eci = [-1.6; 6.8; 0.0001]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β);-sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q)
e = quat2euler(Q)

#state ini
x0 = [(3389.5+125)*1e3; 0.0; 50.0; e; v_eci; 0.0; 0.0; 0.0]
U0 = [50.0;50.0;50.0;0.005;0.005;0.005;1e-3;1e-3;1e-3;1e-3;1e-3;1e-3]
S = [(1e3*3389.5);(1e3*3389.5);(1e3*3389.5); 1.0; 1.0 ;1.0;(1e3*7.00);(1e3*7.00);(1e3*7.00);1.0;1.0;1.0]
S = ones(12)
A0, b0 = ini(x0, U0, S)
t0 = 0.0
tf = 100.0
dtt = 1.0
cond(A0, 2)

Alist, blist, centerlist, XX, T = ellipse_propagation_new(A0, b0, t0, tf, dtt)

Plots.scatter(T, centerlist[1, :]*1e3*3389.5)
Plots.plot!(0.0:0.01:100.0, Z4[1, :])
Plots.scatter(T, centerlist[12, :])
