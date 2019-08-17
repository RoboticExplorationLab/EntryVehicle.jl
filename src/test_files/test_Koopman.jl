#KO : Koopman Operator. Linearize the dynamics at a higher level
using PolyChaos
using DifferentialEquations
using LinearAlgebra
using ODE
using Plots
pyplot()

#generate Fourier-Hermite Polynomial Basis
n = 60 #Number of rec points needed
m = 30 #higher degree of multivariate polynomials (accuracy of the approximation)
op = OrthoPoly("hermite", m, Nrec=n) #if needed Nrec enables to compute more recurrence coefficients
opq = OrthoPolyQ(op, 59)
N = 2 #dimension of the state drives the multivariate number of polynomials
mop = MultiOrthoPoly([opq for i=1:N], m)
P = mop.dim #total number of Polynomials
mop.ind


#compute the tensor products
#t2 = Tensor(2, mop) #needed for G
#t3 = Tensor(3, mop) #needed for A
#t2.get([18,5])

#Use EDMD (Extending Dynamic Mode Decomposition) easier to implement

#First simulate the trajectories (get data)
function pendulum(x_1, dt)
    g = 10.0
    l = 0.1
    ω0 = sqrt(g/l)
    x_2 = zeros(length(x_1))
    x_2[1] = x_1[1]+x_1[2]*dt
    x_2[2] = x_1[2]-(ω0^2)*sin(x_1[1])*dt
    return x_2
end

function propagation(x0, Δt, dt)
    T = 0.0:dt:Δt
    X = zeros(length(x0), length(T))
    x1 = x0
    for i=1:1:length(T)
        x2 = pendulum(x1, dt)
        x1 = x2
        X[:, i] = x2
    end
    X = hcat(x0, X)
    return X
end


#plot(X[1, :], X[2, :])
#plot(X[1, :])

#Construct Matrices:

function matrices(X)
    M = length(X[1, :])
    G = zeros(P, P)
    A = zeros(P, P)
    for k = 1:1:M-1
        Ψ_k1 = [PolyChaos.evaluate(mop.ind[p, :], X[:, k], mop)[1] for p=1:P]
        Ψ_k2 = [PolyChaos.evaluate(mop.ind[p, :], X[:, k+1], mop)[1] for p=1:P]
        G = G + (Ψ_k1)*(Ψ_k1)'
        A = A + Ψ_k1*(Ψ_k2)'
        @show(k)
    end
    G = G/M
    A = A/M
    return A, G
end

x0 = [0.03; 0.1]
X = propagation(x0, 4.0, 0.001)
A, G = matrices(X) #This is the really long part here
K = pinv(G)*A #Koopman Operator
F = eigen(K)
Y = F.vectors
V = F.values

H0 = [PolyChaos.evaluate(mop.ind[p, :], X[:, 1], mop)[1] for p=1:P] #initial value vector
W = (inv(Y))'
Δ = Diagonal(V)
#Δ = Diagonal([log(V[i])/0.01 for i=1:length(V)])
b1 = zeros(P)
b2 = zeros(P)
b1[2] = 1 #because probabilistic Hermite Polynomials
b2[3] = 1
B = hcat(b1, b2)



function KO_evolution(t)
    return (B')*(conj(W)')*exp(Δ*t)*inv(conj(W)')*H0
end

function KO_list()
    t_sim = 0.0:0.01:4
    XX = zeros(2, length(t_sim))
    for i = 1:1:length(t_sim)
        t = t_sim[i]
        XX[:, i] = real(KO_evolution(t))
    end
    return XX
end

x_0 = KO_evolution(0.0)
XX = KO_list()

plot(X[1, :])
plot(XX[1, :])
