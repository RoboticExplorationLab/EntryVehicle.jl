# ##############################################################################
# ################################PCE###########################################
# ##############################################################################

using PolyChaos
using Random
using Distributions

function generate_poly(nrec, d, n, type)
    #type is like "gaussian"
    #nrec number of rec points needed
    #d is the order of the polynomials used
    #n is the dimension of the state vector
    #n = 40 #Number of rec points needed
    #d = 5 #higher degree of multivariate polynomials
    #op = OrthoPoly("gaussian",  d , Nrec=n) #if needed Nrec enables to compute more recurrence coefficients
    #op = GaussOrthoPoly(d)
    opq = GaussOrthoPoly(d; Nrec=nrec) #probabilist Hermite
    mop = MultiOrthoPoly([opq for i=1:n], d)
    return mop
end

#=
P = mop.dim #total number of Polynomials
mop.ind #indices list of for polynomials
showbasis(opq; sym="ξ") #show basis of Polynomials in symbolic language =#

####################################
##############Sampling##############
####################################


function generate_samples_PCE(t0, tf, dt, ξ, Ms, model, n, p)
    t_sim = t0:dt:tf
    samples = zeros(n, length(t_sim), Ms)
    for i=1:1:Ms
        #one column of Z contains the state at a specific time
        t_sim, Z = rk4(model, ξ[:, i], p, dt, [t0; tf])
        samples[:, :, i] = Z
        @show(i)
    end
    return samples
end


####################################
######Coefficients Computation######
####################################

function compute_Phi_matrix(Ms, P, mop, samples, x0, Q0, n, ξ)
    Phi = zeros(Ms, P)
    for i=1:Ms
        for j=1:1:P
            vec = [(ξ[k, i]-x0[k])/sqrt(Q0[k, k]) for k=1:1:n]
            res = PolyChaos.evaluate(mop.ind[j, :], vec, mop)[1]
            Phi[i, j] = res
        end
    end
    return Phi
end

function compute_coeff_PCE_step(samples, t, A, n, P) #t is the time step we want to look at t belongs to
    #function computes PCE coefficients for time step t
    C = zeros(P,n) #contains coefficients for one time step
    for i=1:1:n
        C[:, i] = (A*samples[i, t, :])'
    end
    return C
end

function compute_coeff_PCE_full(samples, Phi, t0, tf, dt, P)
    n, t, Ms = size(samples)
    C = zeros(P, n, t)
    A = pinv(Phi)
    T = t0:dt:tf
    for j = 1:1:length(T)
        c = compute_coeff_PCE_step(samples, j, A, n, P)
        C[:, :, j] = c
        @show(j)
    end
    return C
end

function mean_var_PCE(CC, T, n, P)
    t = length(T)
    m = zeros(n, t)
    var = zeros(n, t)
    for j=1:1:t
        m[:, j] = CC[1, :, j]
        for i = 1:1:n
            var[i, j] = sum(CC[k, i, j]^2 for k=2:1:P)
        end
    end
    return m, var
end

function simulation_PCE(t0, tf, dt, Ms, type, D, d, x0, Q0, model, p)
    #returns the coefficients of the PC expansion for all time steps specified
    #Ms is number of samples (non-intrusive PCE method here)
    #D is the distribution type we want (linked with "type")
    n = length(x0)
    nrec = 100
    mop = generate_poly(nrec, d, n, type)
    DD = MvNormal(x0, Q0/9)
    ξ = rand!(DD, zeros(n, Ms))
    samples = generate_samples_PCE(t0, tf, dt, ξ, Ms, model, n, p) #counting
    P = mop.dim #number of poly in the family
    Phi = compute_Phi_matrix(Ms, P, mop, samples, x0, Q0, n, ξ)
    CC = compute_coeff_PCE_full(samples, Phi, t0, tf, dt, P)
    T = t0:dt:tf
    m, var = mean_var_PCE(CC, T, n, P)
    return T, m, var
end


function sig_pce(var_pce)
    n, t = size(var_pce)
    V = zeros(n, n, t)
    for j = 1:1:t
        V[:, :, j] = Diagonal(var_pce[:, j])
    end
    VV = sig(V)
    return VV
end

#= Example Duffing

t0 = 0.0
tf = 25.0
dt = 0.01
Ms = 2000
type = "gaussian"
D = MvNormal
d = 10
x0 = [0.1;0.1]
Q0 = Diagonal([0.01^2;0.01^2])/9
model = duffing
p =

T_PCE, m_PCE, var_PCE = simulation_PCE(t0, tf, dt, Ms, type, D, d, x0, Q0, model, p)

gr()
VV = [sqrt(var_PCE[1,i]) for i=1:1:length(T_PCE)]
Plots.plot(T_PCE, -3*VV, linestyle = :dashdotdot, color = :green)
Plots.plot!(T_PCE, +3*VV, linestyle = :dashdotdot, color = :green) =#


################################################################################
################## Additional functions 6DOF ###################################
################################################################################

function generate_samples_PCE_6dof(t0, tf, dt, ξ, Ms, model, n, p)
    t_sim = t0:dt:tf
    samples = zeros(n, length(t_sim), Ms)
    for i=1:1:Ms
        #one column of Z contains the state at a specific time
        t_sim, Z = rk4(model, ξ[:, i], p, dt, [t0; tf])
        samples[:, :, i] = point13_12_pce_full(Z)
        @show(i)
    end
    return samples
end

function point12_13_pce(X)
    X13 = zeros(length(X)+1)
    X13[1:3] = X[1:3]
    X13[4:7] = euler2quat(X[4:6])
    X13[8:13] = X[7:12]
    return X13
end

function point13_12_pce(X)
    X12 = zeros(length(X)-1)
    X12[1:3] = X[1:3]
    X12[4:6] = quat2euler(X[4:7])
    X12[7:12] = X[8:13]
    return X12
end

function point13_12_pce_full(Z)
    n, m = size(Z)
    ZZ = zeros(n-1, m)
    for i=1:m
        ZZ[:, i] = point13_12_pce(Z[:, i])
    end
    return ZZ
end

function samples_12_13(ξ)
    n, m = size(ξ)
    ξ13 = zeros(n+1, m)
    for i=1:m
        ξ13[:, i] = point12_13_pce(ξ[:, i])
    end
    return ξ13
end

function simulation_PCE_6dof(t0, tf, dt, Ms, type, D, d, x0_12, Q0_12, model, p)
    #returns the coefficients of the PC expansion for all time steps specified
    #Ms is number of samples (non-intrusive PCE method here)
    #D is the distribution type we want (linked with "type")
    n = length(x0_12)
    nrec = 100
    mop = generate_poly(nrec, d, n, type)
    DD = MvNormal(x0_12, Q0_12/9)
    ξ = rand!(DD, zeros(n, Ms))
    ξ13 = samples_12_13(ξ)
    samples = generate_samples_PCE_6dof(t0, tf, dt, ξ13, Ms, model, n, p) #counting
    P = mop.dim #number of poly in the family
    Phi = compute_Phi_matrix(Ms, P, mop, samples, x0_12, Q0_12, n, ξ)
    CC = compute_coeff_PCE_full(samples, Phi, t0, tf, dt, P)
    T = t0:dt:tf
    m, var = mean_var_PCE(CC, T, n, P)
    return T, m, var
end
