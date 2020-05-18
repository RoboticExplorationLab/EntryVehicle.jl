using LinearAlgebra
using SumOfSquares
using DynamicPolynomials
using Mosek, MosekTools
using JuMP, CSDP
using ForwardDiff
using Plots


################################################################################
############ First Implementation (Tests and comparison to Matlab) #############
################################################################################


function optim_sos_multipliers(var, V, f, h, l1, l2, d1, d2)
    model = SOSModel(with_optimizer(CSDP.Optimizer))
    @variable(model, s1, Poly(monomials(X, 0:d1)))
    @variable(model, s2, Poly(monomials(X, 0:d2)))
    dVdt = differentiate(V, X)'*f
    D = -(dVdt+l2)+s1*(V-1.0)
    E = (h-var)*s2+(1.0-V)
    c1 = @constraint(model, s1 in SOSCone())
    c2 = @constraint(model, s2 in SOSCone())
    c3 = @constraint(model, D in SOSCone())
    c4 = @constraint(model, E in SOSCone())
    optimize!(model)
    return model, c1, c2
end

# test optim_sos_multipliers

model_res, c1, c2 = optim_sos_multipliers(1.2, V, f, h, l1, l2, d1, d2)
termination_status(model_res)


# Bisection method on selecting the t value while solving for polynomials
# s1 and s2 (sos multipliers)

function loop_multipliers_solver(m, V, f, h, l1, l2, d1, d2)
    upper = m+10.0
    lower = m-1.0
    t = 0.0
    s1 = 0.0
    s2 = 0.0
    while abs(upper-lower)>1e-3
        t = (upper+lower)/2
        @show(t)
        model_res, c1, c2 = optim_sos_multipliers(t, V, f, h, l1, l2, d1, d2)
        status = termination_status(model_res)
        @show(status)
        if status == MOI.OPTIMAL || status == MOI.ALMOST_OPTIMAL
            lower = t
            G1 = gram_matrix(c1)
            G2 = gram_matrix(c2)
            s1 = G1.x'*G1.Q*G1.x
            s2 = G2.x'*G2.Q*G2.x
        else
            upper = t
        end
    end
    m = t
    return m, s1, s2
end

function optim_sos_volume(s1, s2, f, h, l1, l2, d)
    model = SOSModel(with_optimizer(CSDP.Optimizer))
    @variable(model, V, Poly(monomials(X, 0:d)))
    @variable(model, β)
    @objective(model, Max, β)
    dVdt = differentiate(V, X)'*f
    D = -(dVdt+l2)+s1*(V-1.0)
    E = (h-β)*s2+(1.0-V)
    c1 = @constraint(model, V-l1 in SOSCone())
    c2 = @constraint(model, D in SOSCone())
    c3 = @constraint(model, E in SOSCone())
    optimize!(model)
    return model, c1
end


function entire_solver(m, k, V, f, h, l1, l2, d1, d2)
    for i=1:k
        m, s1, s2 = loop_multipliers_solver(m, V, f, h, l1, l2, d1, d2)
        @show(m)
        model_resu, c1 = optim_sos_volume(s1, s2, f, h, l1, l2, d)
        m = objective_value(model_resu)
        G = gram_matrix(c1)
        V = G.x'*G.Q*G.x
        @show(m)
    end
    return m, V
end

m, V = entire_solver(0.0, 8, V, f, h, l1, l2, d1, d2)

# tests solvers step by step

@polyvar x1 x2
ϵ = 1e-3
d = 6
d1 = 4
d2 = 2
X = [x1; x2]
f = [-x2; x1+((x1^2)-1)*x2]   # Dynamics of Van Der Pol
A = [0.0 1.0; -1.0 -1.0]   # Transpose of linearized dynamics
Q = [1.0 0.0; 0.0 1.0]
S = lyap(A, Q)
V = (0.65217*x1^2-0.43478*x1*x2+0.43478*x2^2)  # Initial guess V
h = X'*X
l1 = 1e-6*h
l2 = 1e-6*h
m = 0.0
m, s1, s2 = loop_multipliers_solver(m, V, f, h, l1, l2, d1, d2)
model_resu, c1 = optim_sos_volume(s1, s2, f, h, l1, l2, d)


m, V = entire_solver(0.0, 10, V, f, h, l1, l2, d1, d2)

# Plot final level set of function V

x1 = -3.0:0.01:3.0
x2 = -3.0:0.01:3.0
contour(x1,x2,(x1,x2)->V(x1,x2), levels=[1])


################################################################################
########################## Nice implementation  ################################
################################################################################

mutable struct RoaEstimationSolver{T}
    # Region of Attraction estimator
    V::Polynomial{true,T}  # Lyapunov function to be iterated
    s1::Polynomial{true,T}  # sos multiplier
    s2::Polynomial{true,T}  # sos multiplier
    f::Vector{Polynomial{true,T}}  # dynamics vector
    β::T
    X::Vector{PolyVar{true}}
end

struct SolverOptions{T}
    iter::Int  # Number of passes (improvement of solution)
    d1::Int  # Degree of sos multiplier s1
    d2::Int  # Degree of sos multiplier s2
    d::Int  # Degree of Lyapuniv function to find
    ϵ::T  # Sensitivity Bisection Method
    h::Polynomial{true,T}  # Shape factor function (can be modified)
    l1::Polynomial{true,T}  # Used for sos conditioning
    l2::Polynomial{true,T}  # Used for sos conditioning
end

# Not used so far
function Base.copy(s::RoaEstimationSolver{T}) where T
    fnames = fieldnames(typeof(s))
    args = [deepcopy(getfield(s,fname)) for fname in fnames] #deep copy here
    RoaEstimationSolver{T}(args...)
end

function Base.copy(s::SolverOptions{T}) where T
    fnames = fieldnames(typeof(s))
    args = [deepcopy(getfield(s,fname)) for fname in fnames]
    SolverOptions{T}(args...)
end

function optim_sos_multipliers(var, solver::RoaEstimationSolver, opts::SolverOptions)
    l1, l2 = opts.l1, opts.l2
    d1, d2 = opts.d1, opts.d2
    h = opts.h
    V = solver.V
    f = solver.f
    model = SOSModel(with_optimizer(CSDP.Optimizer))
    @variable(model, s1, Poly(monomials(X, 0:d1)))
    @variable(model, s2, Poly(monomials(X, 0:d2)))
    dVdt = differentiate(V, X)'*f
    D = -(dVdt+l2)+s1*(V-1.0)
    E = (h-var)*s2+(1.0-V)
    c1 = @constraint(model, s1 in SOSCone())
    c2 = @constraint(model, s2 in SOSCone())
    c3 = @constraint(model, D in SOSCone())
    c4 = @constraint(model, E in SOSCone())
    optimize!(model)
    return model, c1, c2
end

function step1(solver::RoaEstimationSolver, opts::SolverOptions)
    # Step1 performs the bisection method on the bilinear problem to find the max
    # value of Beta that satisfies the sos problem while finding the sos multipliers.
    # In Step1 the value of th function V is fixed (Lyapunov function is fixed)
    β = solver.β
    upper = β+10.0 # upper bound start
    lower = β-1.0 # lower bound start
    t = 0.0  # temp variable for bisection method
    s1 = solver.s1
    s2 = solver.s2
    while abs(upper-lower)>opts.ϵ
        t = (upper+lower)/2
        model_res, c1, c2 = optim_sos_multipliers(t, solver, opts)
        status = termination_status(model_res)
        if status == MOI.OPTIMAL || status == MOI.ALMOST_OPTIMAL
            lower = t
            G1 = gram_matrix(c1)
            G2 = gram_matrix(c2)
            s1 = G1.x'*G1.Q*G1.x   # Find a way to clean coefficients
            s2 = G2.x'*G2.Q*G2.x   # Find a way to clean coefficients
        else
            upper = t
        end
    end
    β = t
    return β, s1, s2
end

function step2(solver::RoaEstimationSolver, opts::SolverOptions)
    # Step2 optimizes the function V as well as the volume of the estimated region
    s1, s2 = solver.s1, solver.s2
    f = solver.f
    h = opts.h
    l1, l2 = opts.l1, opts.l2
    d = opts.d
    model = SOSModel(with_optimizer(CSDP.Optimizer))
    @variable(model, V, Poly(monomials(X, 0:d)))
    @variable(model, γ)
    @objective(model, Max, γ)
    dVdt = differentiate(V, X)'*f
    D = -(dVdt+l2)+s1*(V-1.0)
    E = (h-γ)*s2+(1.0-V)
    c1 = @constraint(model, V-l1 in SOSCone())
    c2 = @constraint(model, D in SOSCone())
    c3 = @constraint(model, E in SOSCone())
    optimize!(model)
    return model, c1
end

function roa_solve!(solver::RoaEstimationSolver, opts::SolverOptions)
    k = opts.iter
    for i=1:k
        m, m1, m2 = step1(solver, opts)
        solver.s1 = m1
        solver.s2 = m2
        solver.β = m
        model_res, c1 = step2(solver, opts)
        solver.β = objective_value(model_res)
        G = gram_matrix(c1)
        solver.V = (G.x'*G.Q*G.x)+opts.l1
    end
end


function plot_roa!(solver::RoaEstimationSolver)
    V = solver.V
    plt = plot()
    y1 = -3.0:0.01:3.0
    y2 = -3.0:0.01:3.0
    contour!(y1,y2,(y1,y2)->V(y1,y2), levels=[1],
                linewidth=2.0,
                ylabel="X2",
                xlabel="X1",
                title="Van Der Pol ROA Estimation",
                colorbar=:none)
    display(plt)
    return nothing
end


# Test structured solver


#Solver definition and initialization
@polyvar x1 x2
ϵ = 1e-3
d = 6
d1 = 4
d2 = 2
X = [x1; x2]
f = [-1.0*x2; 1.0*x1+((x1^2)-1)*x2]   # Dynamics of Van Der Pol
#A = [0.0 1.0; -1.0 -1.0]   # Transpose of linearized dynamics
#Q = [1.0 0.0; 0.0 1.0]
#S = lyap(A, Q)
V = (0.65217*x1^2-0.43478*x1*x2+0.43478*x2^2)  # Initial guess V
h = 1.0*X'*X
l1 = 1e-6*h
l2 = 1e-6*h
β = 0.0
iter = 8

solver = RoaEstimationSolver(V, x1+0.1*x2, x1+0.1*x2, f, β, X)
solver_opts = SolverOptions(iter, d1, d2, d, ϵ, h, l1, l2)

roa_solve!(solver, solver_opts)
plot_roa!(solver)
