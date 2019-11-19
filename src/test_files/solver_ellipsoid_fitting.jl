using LinearAlgebra
using Random
using Distributions

#Solve Primal
function gradient(H, X)
    #X is the concatenation of the points augmented by 1
    return -inv(H)
end

function second_derivative(H, E1, E2)
    #E1 and E2 are symmetric
    H_inv = inv(H)
    return tr(H_inv*E1*H_inv*E2)
end

#Solve Dual
function H(u, X)
    U = Diagonal(u)
    return X*U*X'
end

function gradient(u, X, m)
    #X is the concatenation of the points augmented by 1
    return [X[:, i]'*H(u, X)*X[:, i] for i=1:1:m]
end

function hessian(u, X, m)
    HESS = zeros(m, m)
    for i = 1:1:m
        for j = 1:1:m
            HESS[i, j] = -(X[:, i]'*H(u, X)*X[:, j])^2
        end
    end
    return HESS
end

function newtons_method(u, ϵ, k_max, X, m)
    k, Δ = 1, fill(Inf, length(u))
    while norm(Δ) > ϵ && k <= k_max
        Δ = inv(hessian(u, X, m))*gradient(u, X, m)
        u -= Δ
        k += 1
        @show(k)
    end
    return u
end


D1 = Uniform(-3, 10)
D2 = Uniform(1, 5)
x1 = zeros(1, 10)
x2 = zeros(1, 10)
rand!(D1, x1)
rand!(D2, x2)
X = vcat(x1, x2, ones(10)')

n, m = size(X);

#init
ϵ = 10^(-9)
k_max = 10000
u = [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0; 10.0]
u = newtons_method(u, ϵ, k_max, X, m)
