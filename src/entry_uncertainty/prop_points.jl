
function prop_points(X, t, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        Xnew[:, i] = dyna_dis(t, dt, X[:, i], u, w)
    end
    return Xnew
end

function prop_points_continuous(X, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, XX = ode78(dyna_coeffoff, X[:, i], 0:0.0001:dt, points=:specified)
        Xnew[:, i] = XX[end]
    end
    return Xnew
end

#this version of points propagation using DiffEq and aero coefficients computed
#offilne.

function prop_points_last(X, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        t_sim, Z = integration2(dyna_coeffoff_inplace, X[:, i], dt)
        Xnew[:, i] = Z[:, end]
    end
    return Xnew
end
