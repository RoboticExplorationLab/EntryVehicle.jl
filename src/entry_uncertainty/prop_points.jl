
function prop_points(X, t, dt, u, w)
    m = length(X[1, :])
    Xnew = zeros(size(X))
    for i=1:1:m
        Xnew[:, i] = dyna_dis(t, dt, X[:, i], u, w)
    end
    return Xnew
end
