#Integrator
using LinearAlgebra
using DifferentialEquations
using Plots
pyplot()

function integration_taylor(I, w_0, t_end, dt)

    C1 = (I[2,2]-I[3,3])/I[1,1]
    C2 = (I[3,3]-I[1,1])/I[2,2]
    C3 = (I[1,1]-I[2,2])/I[3,3]

    A = zeros(12, 12)
    A[1,6] = C1
    A[2,5] = C2
    A[3,4] = C3
    A[4,8] = C2
    A[4,11] = C1
    A[5,7] = C3
    A[5,12] = C1
    A[6,9] = C3
    A[6,10] = C2

    Phi = exp(A*dt)
    T = 0.0:dt:t_end
    W = zeros(3, length(T))
    W[:, 1] = w_0
    for i=1:1:length(T)-1
        W[:, i+1] = Phi[1:3, :]*basis(W[:, i])
    end
    return W
end

function basis(w)
    x = [w;
     w[1]*w[2];
     w[1]*w[3];
     w[2]*w[3];
     w[1]*w[1]*w[2];
     w[1]*w[1]*w[3];
     w[1]*w[2]*w[2];
     w[1]*w[3]*w[3];
     w[2]*w[2]*w[3];
     w[2]*w[3]*w[3]];
     return x
end

#test and comparison with Autodiff

w_0 = [0.1 0.1 0.2]
I = Diagonal([100.0; 200.0; 500.0])
t_end = 30.0
dt = 0.1
T = 0.0:dt:t_end

W = integration_taylor(I, w_0, t_end, dt)
Plots.plot(T, W[1, :])

function euler!(dw, w, p, t)
    I1 = p[1]
    I2 = p[2]
    I3 = p[3]
    C1 = (I2-I3)/I1
    C2 = (I3-I1)/I2
    C3 = (I1-I2)/I3
    dw[1] = C1*w[2]*w[3]
    dw[2] = C2*w[1]*w[3]
    dw[3] = C3*w[1]*w[2]
end

prob = ODEProblem(euler!,w_0,(0.0,t_end), [100.0; 200.0; 500.0])
sol = solve(prob, RK4(), reltol=1e-10,abstol=1e-10)
plot!(sol, vars=(1))


function Koopman(I)
    C1 = (I[2,2]-I[3,3])/I[1,1]
    C2 = (I[3,3]-I[1,1])/I[2,2]
    C3 = (I[1,1]-I[2,2])/I[3,3]

    A = zeros(12, 12)
    A[1,6] = C1
    A[2,5] = C2
    A[3,4] = C3
    A[4,8] = C2
    A[4,11] = C1
    A[5,7] = C3
    A[5,12] = C1
    A[6,9] = C3
    A[6,10] = C2
    return(A)
end

K = Koopman(I)
F = eigen(K)
V = F.values
W = F.vectors
