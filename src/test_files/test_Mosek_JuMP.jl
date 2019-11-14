#Test MOSEK with JuMP
using JuMP
using Mosek
using MosekTools
using Random
using Distributions
using LinearAlgebra
using Convex
using MathOptInterface
using Gurobi

D1 = Uniform(-3, 10)
D2 = Uniform(1, 5)
x1 = zeros(1, 5)
x2 = zeros(1, 5)
rand!(D1, x1)
rand!(D2, x2)
X = vcat(x1, x2)

n, m = size(X);

s = MathOptInterface.LogDetConeTriangle(n)

model = Model(with_optimizer(Mosek.Optimizer, INTPNT_CO_TOL_DFEAS=1e-10))
@variable(model, A[1:n, 1:n])
@variable(model, b[1:n])
#@constraint(model, con[i = 1:m], norm(A*X[:, i]+b)<=1) #no longer supported here
#@NLobjective(model, Max, func(A))
@objective(model, Max, tr(A))
#@constraint(model, A[1, :] in s)
@constraint(model, con[i = 1:m], [1.0; A*X[:, i]+b] in SecondOrderCone())
@constraint(model, A[1, :] in s)

JuMP.optimize!(model)
JuMP.termination_status(model)
A = JuMP.value.(A)
b = JuMP.value.(b)

#plot stuff
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end

ellipse  = A \ B

scatter(X[1,:], X[2,:])
plot!(ellipse[1, :], ellipse[2, :])
