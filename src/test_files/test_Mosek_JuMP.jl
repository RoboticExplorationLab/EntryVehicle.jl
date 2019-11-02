#Test MOSEK with JuMP
using JuMP
using Mosek
using MosekTools
using Random
using Distributions
using LinearAlgebra
using Convex

D1 = Uniform(-3, 10)
D2 = Uniform(1, 5)
x1 = zeros(1, 5)
x2 = zeros(1, 5)
rand!(D1, x1)
rand!(D2, x2)
X = vcat(x1, x2)

n, m = size(X);

model = Model(with_optimizer(Mosek.Optimizer))
@variable(model, A[1:n, 1:n], PSD)
@variable(model, b[1:n])
#@constraint(model, con[i = 1:m], norm(A*X[:, i]+b)<=1) #no longer supported here
@objective(model, Max, logdet(A))
@constraint(model, con[i = 1:m], [1.0; A*X[:, i]+b] in SecondOrderCone())

JuMP.optimize!(model)
JuMP.termination_status(model)
A = value(A)
B = value(B)
