#Linear Covariance Analysis technique

using ForwardDiff

function lincov_prop(x0, Q0, Z, t0, tf, dt, model_lincov)
    #Z is nominal trajectory
    #dt should be same on nominal and desired step here
    #model_lincov is the model modified to integrate forwardDiff
    #E does not make much sense here so far.
    T = t0:dt:tf
    n = length(x0)
    E = zeros(n, length(T))
    QQ = zeros(n, n, length(T))
    x1 = x0
    E[:, 1] = x0
    QQ[:, :, 1] = Q0
    Q1 = Q0
    for i=1:1:length(T)-1
        t = T[i]
        p = [45*pi/180]
        f(u) = model_lincov(t, u, p)
        A = ForwardDiff.jacobian(f, Z[:, i])
        Φ = exp(A*dt)
        x = Φ*x1
        Q = Φ*Q1*Φ'
        E[:, i+1] = x
        QQ[:, :, i+1] = Q
        x1 = x
        Q1 = Q
    end
    return T, QQ
end
