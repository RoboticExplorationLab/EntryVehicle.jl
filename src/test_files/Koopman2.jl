using LinearAlgebra
using DifferentialEquations
using Plots
gr()

#Simple Spectrum LTI system
A = [2.0 0.0; 1.0 3.0]
B = Matrix(adjoint(A))
E = eigen(B)
λ = E.values
W = E.vectors
F = eigen(A)
V = F.vectors

#problem of eigenvalues that are in right order for the matrix A here, that's why
#indices are swapped here after
#just normalized basis of adjoint eigenvectors
W1 = W[:, 1]/(dot(V[:, 2], W[:, 1]))
W2 = W[:, 2]/(dot(V[:, 1], W[:, 2]))

ϕ1(x) = dot(x, W1)
ϕ2(x) = dot(x, W2)

x(x0, t) = exp(λ[1]*t)*ϕ1(x0)*V[:, 2] + exp(λ[2]*t)*ϕ2(x0)*V[:, 1]

x([1.0;1.0], 0.5) #computed using eigenfunctions of KO

#now using regular A matrix
x′(x0, t) = exp(A*t)*x0
x′([1.0;1.0], 0.5)
#--> get the same results, alright that's cool here for x and x′

#More complex non linear system. Whiteboard example

function rk4(f, y_0, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(y_star, t)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(y1, t+dt/2)
        y2 = y_star+k2*dt/2
        k3 = f(y2, t+dt/2)
        y3 = y_star+k3*dt
        k4 = f(y3, t+dt)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y
end

function sys!(du, u, p, t)
    μ = p[1]
    λ = p[2]
    du[1] = μ*u[1]
    du[2] = λ*(u[2]-u[1]^2) #take n is 2 here
end

function sys(u, t)
    du = zeros(2)
    μ = -1.0
    λ = 1.0
    du[1] = μ*u[1]
    du[2] = λ*(u[2]-u[1]^2) #take n is 2 here
    return du
end

p = [-1.0; 1.0]
u0 = [0.1;0.2]
tspan = (0.0, 2.0)

#using DifferentialEquations
prob = ODEProblem(sys!,u0,tspan, p)
sol = DifferentialEquations.solve(prob)
plot(sol, vars=(1, 2))

#using home-made RK4
t_span = [0.0; 2.0]
t_sim, u = rk4(sys, u0, 0.05, t_span)
plot!(u[:, 1], u[:, 2])

#temporal curves
plot(t_sim, u[:, 1])
plot!(t_sim, u[:, 2])


#I now use the solution found by hand for Koopman theory

M = Diagonal([μ; λ; μ*2]) #this is for the fized values of \mu and \lambda
x1(t) = exp(-t)*0.1
H1(t) = exp(t)*(0.2-((0.1)^2)/3)
H2(t) = exp(-2*t)*((0.1)^2)/3

plot(t_sim, [H1(t)+H2(t) for t in t_sim]) #okay we get x2, impressive

X = [0.0 0.0; 1.0 -1.0]
F = eigen(X)
V = F.values
W = F.vectors

W*Diagonal([-1.0; 0.0])*inv(W)

X*[1;1]

inv([0 1; 1 1])
