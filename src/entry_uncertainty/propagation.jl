#propagation file

function propagation(x0, Δt, dt)
    L = 0:dt:Δt
    n = length(x0)
    X = zeros(n, length(L))
    x1 = x0
    u = [0.0]
    for i = 1:1:length(L)
        x2 = dyna_dis_temp(L[i], dt, x1, u, w)
        X[:, i] = x2
        x1 = x2
    end
    return L, X
end
