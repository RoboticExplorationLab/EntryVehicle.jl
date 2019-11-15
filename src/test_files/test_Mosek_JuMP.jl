#Test MOSEK with JuMP
using JuMP
using Mosek
using MosekTools
using Random
using Distributions
using LinearAlgebra
using Convex
using MathOptInterface
using Gurobi #does not support cone exp Primal stuff
using SDPT3
using Plots
using SCS

D1 = Uniform(-3, 10)
D2 = Uniform(1, 5)
x1 = zeros(1, 50)
x2 = zeros(1, 50)
rand!(D1, x1)
rand!(D2, x2)
X = vcat(x1, x2)

n, m = size(X);

s = MathOptInterface.LogDetConeTriangle(n)

model = Model(with_optimizer(Mosek.Optimizer))
@variable(model, A[1:n, 1:n], PSD) #don't forget PSD
@variable(model, b[1:n])
@variable(model , t)
#@constraint(model, con[i = 1:m], norm(A*X[:, i]+b)<=1) #no longer supported here
#@NLobjective(model, Max, func(A))
#@objective(model, Max, tr(A))
#@constraint(model, A[1, :] in s)
@objective(model, Max, t)
@constraint(model, con[i = 1:m], [1.0; A*X[:, i]+b] in SecondOrderCone())
@constraint(model, [t;1;[A[1, 1]; A[1, 2]; A[2, 2]]] in s)

MathOptInterface.TimeLimitSec() = 0.5
JuMP.optimize!(model)
JuMP.termination_status(model)
A = JuMP.value.(A)
b = JuMP.value.(b)
objective_value(model)



#plot stuff
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end

ellipse  = A \ B

scatter(X[1,:], X[2,:])
plot!(ellipse[1, :], ellipse[2, :])

#Comparison with SCS solver
A = Semidefinite(n)
b = Variable(n)
problem = maximize(logdet(A), [norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m])

Convex.solve!(problem, SCSSolver())
problem.optval

A = A.value
b = b.value

#vectorize form of the matrix A

A = [1 2 3; 4 5 6; 7 8 9]
n = 3
V = [A[i, j] for j in 1:1:n for i in 1:1:j]
