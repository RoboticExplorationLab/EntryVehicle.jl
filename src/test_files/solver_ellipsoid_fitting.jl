using LinearAlgebra
using Random
using Distributions
using Plots
using Mosek, Mosek.Ext
using JuMP
using MathProgBase
using MathOptInterface
using MosekTools
gr()

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

#init
ϵ = 10^(-9)
k_max = 10000
u = [1.0; 2.0; 3.0; 4.0; 5.0; 6.0; 7.0; 8.0; 9.0; 10.0]
u = newtons_method(u, ϵ, k_max, X, m)

#Other method from the paper, primal dual method implementation
#reduced optimality conditions

function compute_M_inv_2(u, X)
    n, m = size(X)
    e = ones(m)
    U = Diagonal(u)
    S = 2*(X*U*X'-X*u*u'*X'/(e'*u))
    return S
end


S = compute_M_inv_2(u, X)
S - S'
u = [1.0; 2.0; 3.0; 4.0]

L = Symmetric(S)
L-L'

FF = eigen(L)
W = FF.vectors
D = Diagonal(FF.values)
SS = sqrt(D)
(W*SS*W')^2
cholesky(W*SS*W')

(W*SS*W')-(W*SS*W')' =#

function compute_Σ(u, X)
    n, m = size(X)
    e = ones(m)
    S = compute_M_inv_2(u, X)
    B = (X-X*u*e'/(e'*u))
    Σ = B'*inv(S)*B
    return Σ
end

function jacobian_h(u, X)
    n, m = size(X)
    e = ones(m)
    Σ = compute_Σ(u, X)
    return -2*((Σ/(e'*u))+Σ.*Σ)
end

function compute_h(u, X)
    Σ = compute_Σ(u, X)
    H = [Σ[i, i] for i=1:1:length(Σ[1, :])]
    return H
end

function DRN_direction(u, s, μ, X)
    #u and s must be positive (strictly vectors) and μ positive parameter
    n, m = size(X)
    J = jacobian_h(u, X)
    U = Diagonal(u)
    S = Diagonal(s)
    e = ones(m)
    B = J-inv(U)*S
    r1 = e - s - compute_h(u, X)
    r2 = μ*e-U*s
    Δu = inv(B)*(r1-inv(U)*r2)
    Δs = inv(U)*(r2-S*Δu)
    return Δu, Δs
end

function compute_M(M_inv_2)
    E = inv(M_inv_2)
    @show(E)
    FF = eigen(E)
    W = FF.vectors
    D = Diagonal(FF.values)
    SS = sqrt(D)
    M = W*SS*W'
    return Symmetric(M)
end

function obj(u, X)
    M_inv_2 = compute_M_inv_2(u, X)
    M = compute_M(M_inv_2)
    return -logdet(M)
end

function compute_z(u, X, M)
    n, m = size(X)
    e = ones(m)
    return M*X*u/(e'*u)
end

a=1

function is_non_positive_vector(u)
    C = 0
    m = length(u)
    for i=1:1:m
        if u[i] < 0
            C+=1
        end
    end
    return (C>0)
end

a = 1

function backstep_line_search(u, s, Δu, Δs, X)
    obj_bef = obj(u, X)
    α = 1.0
    u_new = u + α*Δu
    s_new = s + α*Δs
    obj_new = obj(u_new, X)
    @show(obj_bef)
    while #=obj_new >= obj_bef ||=# is_non_positive_vector(u_new) || is_non_positive_vector(s_new)
        α = α/2
        @show(α)
        u_new = u + α*Δu
        s_new = s + α*Δs
        obj_new = obj(u_new, X)
        @show(is_non_positive_vector(u_new))
        @show(is_non_positive_vector(s_new))
        @show(obj_new)
    end
    return α
end

a = 1

function DRN_algo(X)
    n, m = size(X)
    e = ones(m)
    r = 0.99
    α_ini = 1.0
    u_0 = α_ini*ones(m)
    s_0 = α_ini*ones(m)
    u, s = u_0, s_0
    ϵ1 = 10^(-9)
    ϵ2 = 10^(-9)
    i = 0
    while norm(e-compute_h(u, X)-s) > ϵ1 || (u'*s/obj(u, X)) > ϵ2
        @show(i)
        μ = (u'*s)/(10*m)
        Δu, Δs = DRN_direction(u, s, μ, X)
        #@show(u)
        #@show(compute_M(u, X))
        @show(obj(u, X))
        α = backstep_line_search(u, s, Δu, Δs, X)
        α = minimum([r*α; 1.0])
        #α = 0.5
        u = u + α*Δu
        s = s + α*Δs
        @show(is_non_positive_vector(u))
        @show(is_non_positive_vector(s))
        i +=1
    end
    M_inv_2 = compute_M_inv_2(u, X)
    M = compute_M(M_inv_2)
    z = compute_z(u, X, M)
    return M, -z
end

D1 = Uniform(-4, 5)
D2 = Uniform(1.0, 1.0001)
x1 = zeros(1, 15)
x2 = zeros(1, 15)
rand!(D1, x1)
rand!(D2, x2)
X = vcat(x1, x2)
rank([X; ones(15)'])

X = [-4.3 -4.3 2.0 2.0; 1.0 1.0 0.5 0.5]

n, m = size(X);
Plots.scatter(X[1, :], X[2, :])

M, z = DRN_algo(X)

angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) + z[1], sin(angles[i]) + z[2]]
end
ellipse  = M[1:2, 1:2] \ B
Plots.plot!(ellipse[1, :], ellipse[2, :])

function test_DRN(N)
    for i=1:1:N
        D1 = Uniform(-4, 5)
        D2 = Uniform(1.0, 10.0)
        x1 = zeros(1, 15)
        x2 = zeros(1, 15)
        rand!(D1, x1)
        rand!(D2, x2)
        X = vcat(x1, x2)
        M, z = DRN_algo(X)
    end
    return 0.0
end

test_DRN(10000) #okay 4 min for 10000 problem


a=1

function points2ellipse_mosek(X)
    n, m = size(X);
    s = MathOptInterface.LogDetConeTriangle(n)
    model = Model(with_optimizer(Mosek.Optimizer, MSK_DPAR_INTPNT_CO_TOL_DFEAS=10^(-12), MSK_DPAR_INTPNT_CO_TOL_PFEAS=10^(-12), MSK_DPAR_INTPNT_CO_TOL_MU_RED = 10^(-12))) #MSK_DPAR_INTPNT_CO_TOL_INFEAS=10^(-12), MSK_IPAR_INTPNT_MAX_ITERATIONS=1000))
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
    obj = objective_value(model)
    A = JuMP.value.(A)
    b = JuMP.value.(b)
    return A, b, obj
end

A, b = points2ellipse_mosek(X)

angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end
ellipse  = A[1:2, 1:2] \ B
Plots.plot!(ellipse[1, :], ellipse[2, :])

X1 = ellipse2points4(A, b)
Plots.scatter!(X1[1, :], X1[2, :])

A, b = points2ellipse_mosek(X1)
