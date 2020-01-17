using Mosek
using MosekTools
using SCS
using JuMP
using MathOptInterface

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

function prop_points_rk(X, t, dt_rk, p, model, dt_e)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = rk4(model, X[:, i], p, dt_rk, [t, t+dt_e])#integration2(dyna_coeffoff_inplace!, X[:, i], dt)
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

function scale_down_1(X2, S)
    n = length(X2)
    m = 1
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

function uncertainty_propagation(A0, b0, model, t_start, t_end, p, dt_e, dt_rk, scheme, n_scheme)
    T_e = t_start:dt_e:t_end
    n = length(b0)
    blist = zeros(n, length(T_e))
    Alist = zeros(n, n, length(T_e))
    centerlist = zeros(n, length(T_e))
    XX = zeros(n, n_scheme, length(T_e))
    for i=1:1:length(T_e)
        t = T_e[i]
        @show(t)
        X1 = scheme(A0, b0)
        X2 = scale_up(X1, S)
        X3 = prop_points_rk(X2, t, dt_rk, p, model, dt_e)
        X4 = scale_down(X3, S)
        A2, b2 = points2ellipse_mosek(X4)
        blist[:, i] = b0
        Alist[:, :, i] = A0
        centerlist[:, i] = [(-inv(A0)*b0)[j]*S[j] for j=1:n]
        A0 = A2
        b0 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX, T_e
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


function ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p)
    A0, b0 = ini(x0, U0, S)
    n = length(b0)
    Alist, blist, centerlist, XX, T_e = uncertainty_propagation(A0, b0, model, t0, tf, p, dt_e, dt_rk, scheme, n_scheme)
    return Alist, blist, centerlist, XX, T_e
end

function variance_e(Alist, S)
    n, n, t = size(Alist)
    var_e = zeros(size(Alist))
    for i=1:1:t
        var_e[:, :, i] = inv(Alist[:, :, i]*Alist[:, :, i])
    end
    V = sig(var_e, S)
    return V
end


#= Example Duffing
x0 = [0.1; 0.1]
U0 = [0.02; 0.02]
S = [1.0; 1.0]
t0 = 0.0
tf = 25.0
dt_rk = 0.01
dt_e = 0.1
scheme = ellipse2points4
n = length(x0)
n_scheme = 2*n
model = duffing
p = [-1.0; 1.0; 0.2; 0.1; 1.0]

Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p)

Plots.plot(centerlist[1, :], centerlist[2, :]) =#


#= Example Vinhs
δ = 70.0*pi/180.0
r_min = 0.2
r_cone = 1.3
r_G = [0.0; 0.0; -0.189]
A_ref = pi*r_cone^2
table_CD, table_CL = drag_lift_table(δ, r_min, r_cone, r_G)

S = [3389.53*1e3; 1.0; 1.0; 1e3*7.00; 1.0; 1.0]
x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
U0 = [50.0; (0.0017); (0.0017); (1.0); (0.0001); (0.0001)]
t0 = 0.0
tf = 80.0
dt_rk = 0.01
dt_e = 1.0
scheme = ellipse2points
n = length(x0)
n_scheme = 2*n+1
model = vinhs_full_model
p = [45.0*pi/180]

Alist, blist, centerlist, XX, T_e = ellipse_prop(x0, U0, S, t0, tf, dt_rk, dt_e, scheme, n_scheme, model, p)

Plots.plot(T_e, centerlist[6, :])
Plots.plot(T_e, centerlist[2, :]) =#


################################################################################
#################### Additional functions 6DOF #################################
################################################################################

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

function ellipse_propagation(A0, b0, tstart, tend, dtt, x0_i, model, p)
    T = tstart:dtt:tend
    n = length(b0)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(13, 2*n+1, length(T))
    ref = zeros(13, length(T))
    p = [0.0]
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A0, b0, x0_i) #return a set of points in dim 13
        X2 = scale_up_6dof(X1, S)
        X3 = prop_points_rk(X2, t, dt_rk, p, model, dtt)
        x1 = X3[:, end]
        X_12 = points13_12(X3, x1)
        X4 = scale_down_6dof(X_12, S)
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

function scale_up_6dof(X3, S)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[1, i] = X3[1, i]*S[1]
        X4[2, i] = X3[2, i]*S[2]
        X4[3, i] = X3[3, i]*S[3]
        X4[4, i] = X3[4, i]
        X4[5, i] = X3[5, i]
        X4[6, i] = X3[6, i]
        X4[7, i] = X3[7, i]
        X4[8, i] = X3[8, i]*S[7]
        X4[9, i] = X3[9, i]*S[8]
        X4[10, i] = X3[10, i]*S[9]
        X4[11, i] = X3[11, i]*S[10]
        X4[12, i] = X3[12, i]*S[11]
        X4[13, i] = X3[13, i]*S[12]
    end
    return X4
end


function scale_down_6dof(X3, S)
    n, m = size(X3)
    X4 = zeros(n, m)
    for i=1:1:m
        X4[:, i] = [X3[j, i]/S[j] for j=1:length(S)]
    end
    return X4
end

function euler2quat_center(centerlist, ref)
    n, t = size(centerlist)
    center_X13 = zeros(n+1, t)
    for i=1:1:t
        center_X13[:, i] = point12_13(centerlist[:, i], ref[:, i])
    end
    return center_X13
end

function scale_centerlist_up(centerlist, S)
    #centerlist is of dimension 12
    n, t = size(centerlist)
    C = zeros(n, t)
    for j=1:1:t
        C[:, j] = [centerlist[i, j]*S[i] for i=1:length(S)]
    end
    return C
end

function transform_centerlist(centerlist, ref, S)
    C = scale_centerlist_up(centerlist, S)
    center_X13 = euler2quat_center(C, ref)
    return center_X13
end


function scale_down_13(X, S)
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

#additional sampling schemes

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
