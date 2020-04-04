#Extended Dynamic Mode Decomposition
#Simple method for solving the Koopman operator

#The first part of the code does Koopman approximation using the Hermite Pol
#Second part is application of EDMD method to linear and duffing

################################################################################
# Reference: A DATA–DRIVEN APPROXIMATION OF THE KOOPMAN OPERATOR:
# EXTENDING DYNAMIC MODE DECOMPOSITION
################################################################################

using LinearAlgebra
using Plots


# Dynamics #####################################################################
################################################################################

P = [-1.0, 1.0, 0.2, 0.1, 1.0] #Duffing model parameters

function duffing(u,P,t)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = P[1]
    β = P[2]
    δ = P[3]
    γ = P[4]
    ω = P[5]
    du[1] = u[2]
    du[2] =  γ*cos(t)-δ*u[2]-α*u[1]-β*u[1]^3
    return du
end

function duffing_dis(u,P,t,dt)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = P[1]
    β = P[2]
    δ = P[3]
    γ = P[4]
    ω = P[5]
    du[1] = u[1]+u[2]*dt
    du[2] =  u[2]+(γ*cos(t)-δ*u[2]-α*u[1]-β*u[1]^3)*dt
    return du
end

dt = 1e-3
t_end = 100.0
T = 0:dt:t_end
X0 = [0.1, 10]

function prop_dis(T, dt, X0)
    X = zeros(2, length(T))
    X1 = X0
    X[:, 1] = X0
    for i=1:length(T)-1
        X[:, i+1] = duffing_dis(X1,P,T[i],dt)
        X1 = X[:, i+1]
    end
    return X
end

X = prop_dis(T, dt, X0)
Plots.plot(X[1, :], X[2, :])


# Generate Legendre Polynomials ################################################
################################################################################

using PolyChaos

#=
n = 2 #dimension
nrec = 40 #Number of rec points needed
d = 4 #higher degree of multivariate polynomials
opq = GaussOrthoPoly(d; Nrec=nrec) #Legendre
mop = MultiOrthoPoly([opq for i=1:n], d)
P = mop.dim #total number of Polynomials
mop.ind
showbasis(opq; sym="ξ") =#

function generate_hermite(n, p, nrec=40)
    opq = GaussOrthoPoly(p; Nrec=nrec) #Legendre
    mop = MultiOrthoPoly([opq for i=1:n], p)
    return mop
end


function nbr_terms_of_deg(d)
    #returns the ndr of terms of total degree d for a n dimensiona polynomial
    return d+1
end

function nbr_lines_before_deg(d)
    S = 0
    for i=0:1:d-1
        S += nbr_terms_of_deg(i)
    end
    return S
end


function map_index(i, j, p)
    d = i+j #return the degree of the term
    if d == 0
        return 1
    else
        S = nbr_lines_before_deg(d) + (d-i+1)
    end
end

function koopman_hermite_matrix(p, t, P)
    #p is the order at which we truncate the linearization (polynomial order)
    α = P[1]
    β = P[2]
    δ = P[3]
    γ = P[4]
    N = Int64((p+1)*(p+2)/2)
    A = zeros(N, N)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p) #number of the line we are working at
            if k <= N && k >=1
                A[k, k] = -δ*j #ok same
                k1 = map_index(i-1, j+1, p)
                (k1 <= N && k1 >=1) ? A[k, k1] = i : Nothing
                k2 = map_index(i-1, j-1, p)
                (k2 <= N && k2 >=1) ? A[k, k2] = i*j -α*i*j -β*j*3*i^2 : Nothing
                k3 = map_index(i, j-1, p)
                (k3 <= N && k3 >=1) ? A[k, k3] = j*γ*cos(t) : Nothing
                k4 = map_index(i, j-2, p)
                (k4 <= N && k4 >=1) ? A[k, k4] = -δ*j*(j-1) : Nothing
                k5 = map_index(i+1, j-1, p)
                (k5 <= N && k5 >=1) ? A[k, k5] = -α*j-β*j*(3*i+3) : Nothing
                k6 = map_index(i+3, j-1, p)
                (k6 <= N && k6 >=1) ? A[k, k6] = -β*j : Nothing
                k8 = map_index(i-3, j-1, p)
                (k8 <= N && k8 >=1) ? A[k, k8] = -β*j*i*(i-1)*(i-2) : Nothing
            end
        end
    end
    return A
end

p = 3
A = koopman_hermite_matrix(p, 1.0, P)
map_index(4, 0, 4)

#=
function basis(x, p)
    N = Int64((p+1)*(p+2)/2)
    base = zeros(N)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p)
            base[k] = PolyChaos.evaluate(mop.ind[k, :], x, mop)[1]
        end
    end
    return base
end =#

function basis(x, p)
    N = Int64((p+1)*(p+2)/2)
    base = zeros(N)
    n,m = size(mop.ind)
    for k=1:n
        base[k] = PolyChaos.evaluate(mop.ind[k, :], x, mop)[1]
    end
    return base
end

function integration_koopman_hermite(x_ini, t_ini, t_end, dt)
    T = t_ini:dt:t_end
    X = zeros(length(x_ini), length(T))
    X[:, 1] = x_ini
    for q=1:1:length(T)-1
        A = koopman_hermite_matrix(p, T[q], P)
        Phi = exp(A*dt)
        k1 = map_index(1, 0, p)
        k2 = map_index(0, 1, p)
        X[:, q+1] = vcat(Phi[k1, :]', Phi[k2, :]')*basis(X[:, q], p)
    end
    return X
end

# Computation ##################################################################
################################################################################

n = 2 #dimension 3 system
p = 4 #max degree of polynomials
mop = generate_hermite(n, p)
show(mop)
P = [-1.0, 1.0, 0.2, 0.1, 1.0] #Duffing model parameters
t_ini = 0.0
t_end = 100.0
dt = 1e-3
x_ini = [0.1 0.1]
N = Int64((p+1)*(p+2)/2)
X2 = integration_koopman_hermite(x_ini, t_ini, t_end, dt)

T = t_ini:dt:t_end
X = prop_dis(T, dt, x_ini)
Plots.plot(X[1, :], X[2, :])
Plots.plot(X2[1, :], X2[2, :], legend = false)




# EDMD implementation ##########################################################
################################################################################

using Random
using Distributions


function linear_ex(x1)
    J = [0.9 -0.1; 0.0 0.8]
    return J*x1
end


function generate_data(M)
    #generate random data for Koopman learning process
    X = zeros(2, M)
    Y = zeros(2, M)
    for i=1:M
        x1 = randn(1)
        x2 = randn(1)
        x_ini = [x1;x2]
        y = linear_ex(x_ini)
        X[:, i] = x_ini
        Y[:, i] = y
    end
    return X, Y
end

function generate_hermite(n, p, nrec=40)
    opq = GaussOrthoPoly(p; Nrec=nrec) #Legendre
    mop = MultiOrthoPoly([opq for i=1:n], p)
    return mop
end

function basis_hermite(x, p)
    N = Int64((p+1)*(p+2)/2)
    base = zeros(N)
    n,m = size(mop.ind)
    for k=1:n
        base[k] = PolyChaos.evaluate(mop.ind[k, :], x, mop)[1]
    end
    return base
end

function compute_A_G(X, Y, N)
    n,M = size(X)
    A = zeros(N, N)
    G = zeros(N, N)
    for i=1:M
        x_m = X[:, i]
        y_m = Y[:, i]
        Phi_x = basis_hermite(x_m, p)
        Phi_y = basis_hermite(y_m, p)
        A+=Phi_x*Phi_y'
        G+=Phi_x*Phi_x'
    end
    return A, G
end


# The basis of functions depends on the dynamics and the way we generated data !!

#main computation
n = 2
p = 3 # max degree of polynomial basis
M = 10000 #number of data pairs
N = Int64((p+1)*(p+2)/2)
mop = generate_hermite(n, p)
X, Y = generate_data(M)
A, G = compute_A_G(X, Y, N)
K = pinv(G)*A #linear truncated approx of Koopman matrix
F = eigen(K)
V = F.values
W = F.vectors

#one scalar eigenfunction of Koopman is:
function eigenfunctions(x, p, W)
    return W'*basis_hermite(x, p)
end

#Plotting contour plot of eigenfunctions

x = -5.0:0.1:5.0
y = -5.0:0.1:5.0

function plot_eigenfunctions(x, y, W, p)
    # x and y gives the domain on which we plot
    # plot the function with infinite norm equals to 1 on the given domain
    XX = repeat(reshape(x, 1, :), length(y), 1)
    YY = repeat(y, 1, length(x))
    ZZZ = zeros(length(x), length(y), length(W[1, :]))
    P = []
    for i=2:length(W[1,:]) #skip the first eigenfunction (trivial one)
        f(x, y) = begin
            eigenfunctions([x,y], p, W)[i]
        end
        ZZ = map(f, XX, YY)
        nor = maximum(abs.(ZZ))
        ZZZ[:, :, i] = ZZ./nor
    end
    return ZZZ
end

#Plot successively (see how to automate that)
ZZZ = plot_eigenfunctions(x, y, W, p)
i = 9
mop.ind[i, :]
Plots.contourf(x, y, ZZZ[:, :, i])
#note : some of the eigenfunctions might be trivial ones.

# EDMD on Duffing Oscillator ###################################################
################################################################################


function duffing_dis_aut(u,P,dt)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = P[1]
    β = P[2]
    δ = P[3]
    du[1] = u[1]+u[2]*dt
    du[2] =  u[2]+(-δ*u[2]-α*u[1]-β*u[1]^3)*dt
    return du
end

function generate_data(M)
    #generate random data for Koopman learning process
    #might need to change the tirage
    X = zeros(2, M)
    Y = zeros(2, M)
    D = Uniform(-2, 2)
    for i=1:M
        x1 = rand(D)
        x2 = rand(D)
        x_ini = [x1;x2]
        y = duffing_dis_aut(x_ini,P,dt)
        X[:, i] = x_ini
        Y[:, i] = y
    end
    return X, Y
end

function generate_hermite(n, p, nrec=40)
    opq = GaussOrthoPoly(p; Nrec=nrec) #Legendre
    mop = MultiOrthoPoly([opq for i=1:n], p)
    return mop
end

function basis_hermite(x, p)
    N = Int64((p+1)*(p+2)/2)
    base = zeros(N)
    n,m = size(mop.ind)
    for k=1:n
        base[k] = PolyChaos.evaluate(mop.ind[k, :], x, mop)[1]
    end
    return base
end

function compute_A_G(X, Y, N)
    n,M = size(X)
    A = zeros(N, N)
    G = zeros(N, N)
    for i=1:M
        x_m = X[:, i]
        y_m = Y[:, i]
        Phi_x = basis_hermite(x_m, p)
        Phi_y = basis_hermite(y_m, p)
        A+=Phi_x*Phi_y'
        G+=Phi_x*Phi_x'
    end
    return A, G
end


function eigenfunctions(x, p, W)
    return W'*basis_hermite(x, p)
end


function plot_eigenfunctions_duff(x, y, W, p)
    # x and y gives the domain on which we plot
    # plot the function with infinite norm equals to 1 on the given domain
    XX = repeat(reshape(x, 1, :), length(y), 1)
    YY = repeat(y, 1, length(x))
    ZZZ_n = zeros(length(x), length(y), length(W[1, :]))
    ZZZ_p = zeros(length(x), length(y), length(W[1, :]))
    P = []
    for i=2:length(W[1,:]) #skip the first eigenfunction (trivial one)
        f(x, y) = begin
            eigenfunctions([x,y], p, W)[i]
        end
        ZZ = map(f, XX, YY)
        ZZZ_n[:, :, i] = norm.(ZZ)
        ZZZ_p[:, :, i] = angle.(ZZ) #radians
    end
    return ZZZ_n, ZZZ_p
end


# Main computation
dt = 1e-3
n = 2
p = 3 # max degree of polynomial basis
P = [1.0,1.0,0.5] #Duffing parameters matching paper
M = 1000
N = Int64((p+1)*(p+2)/2)
mop = generate_hermite(n, p)
X, Y = generate_data(M)
A, G = compute_A_G(X, Y, N)
K = pinv(G)*A #linear truncated approx of Koopman matrix
F = eigen(K)
V = F.values
W = F.vectors


#Plotting contour plot of eigenfunctions

x = -2.0:0.1:2.0
y = -2.0:0.1:2.0

#Plot successively (see how to automate that)
ZZZ_n, ZZZ_p = plot_eigenfunctions_duff(x, y, W, p)
i = 2
Plots.contourf(x, y, ZZZ_n[:, :, i])
Plots.contourf(x, y, ZZZ_p[:, :, i])

#Actually need to use RBF (using spline method given in paper)
