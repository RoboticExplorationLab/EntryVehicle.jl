using LinearAlgebra
using Plots
using DifferentialEquations
using ODEInterfaceDiffEq
pyplot()
gr()

f(u, p, t) = 4*u*(1-u);
u0 = pi/10;
tspan = (0.0, 100.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, AutoTsit5(Rosenbrock23()), reltol=1e-10,abstol=1e-10)

plot(sol.t, sol.u)
xlabel!("time [s]")
ylabel!("position [m]")
title!("Time Evolution")
#Does not produce what I would expect. Might be the ODE solver of something else
#OKAY solved. This is due to the difference between discrete and continuous systems...
#Discrete dynamic system here.

function discrete_dyna(t, x_1)
    x_2 = 4*x_1*(1-x_1)
    return x_2
end

function prop(x_0)
    time = 1:1:100
    X = zeros(length(time)+1)
    X[1] = x_0
    x_1 = x_0
    for i = 1:length(time)
        X[i+1] = discrete_dyna(time[i], x_1) #might need to tweak time here I think if dyna depends on t
        x_1 = discrete_dyna(time[i], x_1)
    end
    t_sim = 0:1:100
    return t_sim, X
end

#See the sensitivity to initial conditions on this example
t_sim, X = prop(pi/10)
plot(t_sim, X)

t_sim2, X2 = prop(pi/10+0.001)
plot(t_sim2, X2)

##############################################
############ KOOPMAN DECOMPOSITION ###########
##############################################
#test Koopman decomposition
using PolyChaos
using LinearAlgebra
using Plots
pyplot()

#get hermite polynomials in 1D
n = 60 #Number of rec points needed
m = 10 #higher degree of multivariate polynomials (accuracy of the approximation)
op = OrthoPoly("gaussian", m, Nrec=n) #if needed Nrec enables to compute more recurrence coefficients
#Gaussian is the regular probabilistic hermite polynomials here.
opq = OrthoPolyQ(op, 59)
N = 1 #dimension of the state drives the multivariate number of polynomials
mop = MultiOrthoPoly([opq for i=1:N], m)
P = mop.dim #total number of Polynomials
mop.ind

t2 = Tensor(2, mop)
t3 = Tensor(3, mop)

#okay nice here goes from 0 to 10 (makes sense)
t2.get([0, 0]) #gives the norm (basically it is the factorial here)

#so by and large, this is my dynamics here in fact
F = zeros(m+1);
F[1:3] = [1.0, -1.0, 1.0];

function H(x)
    return [PolyChaos.evaluate(mop.ind[p, :], x, opq)[1]/t2.get([p-1, p-1]) for p=1:P] #not sure about the normalization here
end

function K_finite()
    G = zeros(m+1, m+1) #but n and m should be different here.
    A = zeros(m+1, m+1)
    for i=1:m+1
        G[i, i] = t2.get([i-1, i-1]) #important shift here
        for j=1:m+1
            A[i, j] = (j-1)*sum(F[k]*t3.get([i-1, k-1, j-1]) for k=1:m+1)
        end
    end
    return inv(G)*A
end

#construct the finite approximation of the Koopman operator or similarly the projection of the operator on the 11-dim space
K = K_finite()
#extract the eigenvalues and vectors in this finite space (coordinates for the functions basis)
V = eigvals(K)
Σ = eigvecs(K)
Σ*Diagonal(V)*inv(Σ) #Okay this gives back K so alright
W = inv(Σ)'
Δ = Diagonal(V)
Q = inv(eigvecs(K)) #supposedly left eigenvectors
W = Q
W'*Σ

for i=1:length(Q[:, 1]) #matrix of left eigenvectors. Properly scaled?
    N = Q[:, i]'*Σ[:, i]
    Q[:, i] = Q[:, i]/N
end
Q[:, 5]'*Σ[:, 5]


#coordinates of the identity observable on the hermite basis
B = zeros(m+1);
B[2] = 1.0;

function u_Koopman(t, u0)
    K = K_finite()
    V = eigvals(K)
    Σ = eigvecs(K) #matrix of right eigenvectors
    Q = inv(eigvecs(K))
    for i=1:length(Q[:, 1]) #matrix of left eigenvectors. Properly scaled?
        N = Q[:, i]'*Σ[:, i]
        Q[:, i] = Q[:, i]/N
    end
    Δ = Diagonal(V)
    B = zeros(m+1);
    B[2] = 1.0;
    u = B'*Q'*exp(Δ*t)*inv(Q')*H(u0)
    return u
end

H(5.0) #might be that this function should be normalized
u_Koopman(0.5, 0.1) #big values. Something wrong somewhere


function exact_solution(t, u0)
    return u0/(u0+(1-u0)*exp(t))
end

#Exact solution plotting
T = 0.0:0.1:10
plot(T, [exact_solution(t, 0.1) for t in T])

#Koopman Approximation plotting
plot(T, [u_Koopman(t, 0.1) for t in T])

#dont think I need to scale the H thing otherwise it does not make sense at t=0
#but might be a scaling issue when defining Q
#might also be that the Hermite basis is not adapted
#also the dimensions might no be enough

#Discretized version of the dynamics
function discrete(t, x1)
    x2 = -x1+x1^2
    return x2
end

function prop_K(x_0)
    time = 1:1:100
    X = zeros(length(time)+1)
    X[1] = x_0
    x_1 = x_0
    for i = 1:length(time)
        X[i+1] = discrete(time[i], x_1) #might need to tweak time here I think if dyna depends on t
        x_1 = discrete(time[i], x_1)
    end
    t_sim = 0:1:100
    return t_sim, X
end

t_sim, X = prop_K(0.5)

#Plot discrete version of dynamics (not discretized but discrete!)
plot(t_sim, X)


#okay so try integrator on the system du/dt = -u+u^2

function basis(u)
    x = [u;u^2;u^3;u^4;u^5;u^6;u^7;u^8;u^9;u^10]
    return x
end

function integrator(u0, tend, dt)
    A = zeros(10, 10)
    for i=1:1:9
        A[i,i] = -i
        A[i,i+1] = i
    end
    A[10,10] = -10
    Phi = exp(A*dt)
    T = 0.0:dt:t_end
    U = zeros(length(T))
    U[1] = u0
    for i=1:1:length(T)-1
        U[i+1] = Phi[1, :]'*basis(U[i])
    end
    return T, U
end

function exact_solution(t, u0)
    return u0/(u0+(1-u0)*exp(t))
end

u0 = 0.1
t_end = 30.0
dt = 0.1
t_sim, U = integrator(u0, t_end, dt)

plot(t_sim, U)
plot!(t_sim, [exact_solution(t, u0) for t in t_sim])

V = [exact_solution(t, u0) for t in t_sim]
W = zeros(length(t_sim))
for i=1:1:length(U)
    W[i] = U[i]-V[i]
end
W

plot(t_sim, W)

#test to diagonalize matrix A here
A = zeros(10, 10)
for i=1:1:9
    A[i,i] = -i
    A[i,i+1] = i
end
A[10,10] = -10

F = eigen(A)
V = F.values
W = F.vectors

W*Diagonal(V)*inv(W) #okay that gives A back (almost)

inv(W)
