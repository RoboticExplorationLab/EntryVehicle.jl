# Test Sum Of Squares julia

using LinearAlgebra
using SumOfSquares
using DynamicPolynomials
using Mosek
using MosekTools
using JuMP
using CSDP
using MathProgBase
using MathOptInterface

# first trial using example from DOC ###########################################

# Create symbolic variables (not JuMP decision variables)
@polyvar x1 x2

# Create a Sum of Squares JuMP model with the Mosek solver
model = SOSModel(with_optimizer(Mosek.Optimizer))

# Create a JuMP decision variable for the lower bound
@variable(model, γ)

# f(x) is the Goldstein-Price function
f1 = x1+x2+1
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2
f3 = 2*x1-3*x2
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2

f = (1+f1^2*f2)*(30+f3^2*f4)

# Constraints f(x) - γ to be sum of squares
@constraint(model, f >= γ)

@objective(model, Max, γ)

optimize!(model)

# The lower bound found is 3
println(objective_value(model))
################################################################################

# Example Risk Aware and Robust Nonlinear Planning MIT  ########################

@polyvar x1
p = x1^4+4*x1^3+6*x1^2+4*x1+5
model = SOSModel(with_optimizer(Mosek.Optimizer))
@variable(model, γ)
cref = @constraint(model, p >= γ)
@constraint(model, γ>= 0.0)
@objective(model, Min, γ)
optimize!(model)
println(objective_value(model))

# Attempt to determine if this is SOS more rigorously ##########################
model2 = SOSModel(with_optimizer(Mosek.Optimizer)) #feasibility pb by default.
@constraint(model2, p>=0.0)
@constraint(model2, p in SOSCone()) #okay does the same
optimize!(model2)
println(objective_value(model2)) # 0.0 ofc because feasibility problem
println(primal_status(model2)) #this gives "FEASIBILITY POINT" or "NO SOLUTION"
println(dual_status(model2))
println(termination_status(model2))

termination_status(model2) == MOI.OPTIMAL

# Attempt with constraint on the region ########################################
@polyvar x y
model3 = SOSModel(with_optimizer(Mosek.Optimizer))
q = x^3-4*x^2+2*x*y-y^2+y^3
S = @set x>=0.0 && y>=0.0 && x+y-1.0>=0.0
@constraint(model3, q>=0.0, domain=S) # Actually it searches if it is SOS,
#does not check nonnegativity but in this case it is equivalent I believe
optimize!(model3)
println(primal_status(model3))
println(dual_status(model3))

all_constraints(model3, VariableRef, MOI.GreaterThan{Float64})
sos_decomposition(p in SOSCone())
constraint_ref_with_index(model2,2)



@constraint(model, p in SOSCone())
# Comparison with MATLAB code MIT class on robust planning

@polyvar x1 x2
p = x1^2 -x1*x2^2 +x2^4+1
model = SOSModel(with_optimizer(CSDP.Optimizer))
X = monomials([x1, x2], 0:2)
n = length(X)
cref = @constraint(model, p in SOSCone())
gram_matrix(@constraint(model, p>=0.0))
@variable(model, q, SOSPoly(X))
@constraint(model, p == q)
optimize!(model)
Q = value(q)[:, :]
X'*Q*X == p
cholesky(Q)
println(objective_value(model))
println(p)
gram_matrix(p)


Q = gram_matrix(cref)
Q.Q[:, :]

Q.x'*Q.Q*Q.x
certificate_monomials(cref)
sos_decomposition(cref, 1e-4)

gram_matrix(cref).x


status(model2)

#Lyapunov Function Example 1 MIT class #########################################
################################################################################
using SumOfSquares
using DynamicPolynomials
using Mosek
using MosekTools
using CSDP
using Plots

@polyvar x1 x2
f=[-x1+(1+x1)*x2;-(1+x1)*x1]  #dynamics equation
x = [x1;x2]
model = SOSModel(with_optimizer(CSDP.Optimizer))

d = 4 #degree of the Lyapunov function that we seek here
@variable(model, V, Poly(monomials(x, 0:d)))
dVdt = differentiate(V, x)'*f
c1 = @constraint(model, V in SOSCone())
c2 = @constraint(model, -dVdt in SOSCone())
c3 = @constraint(model, V[length(V)]==0)
optimize!(model)

lyap = value(V)

y1 = -1.0:0.01:1.0
y2 = -1.0:0.01:1.0
surface(y1,y2,(y1,y2)->lyap(y1,y2),
            linewidth=1.0,
            ylabel="X2",
            xlabel="X1",
            title="Quartic lyapunov function", colorbar=:none)
