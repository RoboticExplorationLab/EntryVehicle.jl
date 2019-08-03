
function prop_points(X, t, dt, u, w)
    n = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:n
        Xnew[:, i] = dyna_dis(t, dt, X[:, i], u, w)
    end
    return Xnew
end
