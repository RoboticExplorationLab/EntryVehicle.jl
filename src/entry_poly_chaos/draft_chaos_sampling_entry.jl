#Using Polynomial Chaos for entry uncertainty without control
#Sampling Approach (Regression)
using LinearAlgebra
using PolyChaos
using ODE
using DifferentialEquations
using Statistics
using Distributions
using Random
using Plots
using KernelDensity

pyplot()

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model.jl")
include("propagation.jl")
include("animate_sampling.jl")

####################################
###########Generate Basis###########
####################################
#assume Gaussian disturbances ==> Probabilistic Hermite Polynomials

n = 40 #Number of rec points needed
d = 5 #higher degree of multivariate polynomials
op = OrthoPoly("gaussian",  d , Nrec=n) #if needed Nrec enables to compute more recurrence coefficients
#op = GaussOrthoPoly(d)
opq = GaussOrthoPoly(d; Nrec=n) #probabilist Hermite
N = 6 #number of random inputs
mop = MultiOrthoPoly([opq for i=1:N], d)
P = mop.dim #total number of Polynomials
mop.ind
showbasis(opq; sym="ξ")



####################################
##############Sampling##############
####################################

Ms = 2000 #number of sampled trajectories
Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
q = mat2quat(M)
q = qconj(q)
x0 = [(3389.5+125)/Re, 0.0, 0.0, q[1], q[2], q[3], q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Q0 = Diagonal([1.0/Re; 1.0/Re; 1.0/Re; 0.0000001; 0.000001; 0.000000001; 0.0000001; 0.001; 0.001; 0.001; 0.001; 0.001; 0.001])
Q0 = Matrix(Q0) #Covariance Matrix (no correlation)

D = MvNormal(x0, Q0)
ξ = zeros(length(x0), Ms)
rand!(D, ξ) #one column is one state sample in x

#integration parameters
t0 = 0.0
dt = 1.0 #stored every unit (no implication on solver)
tf = 50.0
t_sim = t0:dt:tf
w = [0.0158*10^9; 0.0; 0.0; 0.0] #can be varied later

δ = 70*pi/180
r_cone = 1.3
r_G = [0.2; 0.0; 0.3]
table_CF, table_Cτ = table_aero(δ, r_cone, r_G) #offline coefficients computation

samples = zeros(13, length(t_sim), Ms)

for i=1:1:Ms
    function dyn(t, x)
        return dyna_coeffoff(t, x, [0.0], w) #w can be varied later
    end
    Z = propagation(dyn, ξ[:, i], t_sim) #one column of Z contains the state at a specific time
    samples[:, :, i] = Z
    @show(i)
end

####################################
######Coefficients Computation######
####################################

function compute_Phi_matrix(Ms, P, mop, samples, x0, Q0)
    Phi = zeros(Ms, P)
    for i=1:Ms
        for j=1:1:P
            vec = [(ξ[k, i]-x0[k])/sqrt(Q0[k, k]) for k=1:1:13]
            res = PolyChaos.evaluate(mop.ind[j, :], vec, mop)[1]
            Phi[i, j] = res
        end
    end
    return Phi
end

function compute_coeff_PCE(samples, t, Phi) #t is the time step we want to look at t belongs to
    #function computes PCE coefficients for time step t
    C = zeros(P, 13) #contains coefficients for one time step
    for i=1:1:13
        C[:, i] = pinv(Phi)*samples[i, t, :]
    end
    return C
end

Phi = compute_Phi_matrix(Ms, P, mop, samples, x0, Q0)
C = compute_coeff_PCE(samples, 50, Phi) #each column contains the coeff PCE of each state variable at time t


#reconstruct density distributions (detail for x position)
Nr = 1000 #samples for reconstruction
Ξ = zeros(length(x0), Nr)
rand!(D, Ξ)
function compute_u()
    ũ_x = zeros(Nr, 1)
    for p=1:1:P
        ũ_x = ũ_x + [C[p, 1]*PolyChaos.evaluate(mop.ind[p, :], Ξ[:, i], mop)[1] for i=1:Nr]
    end
    return ũ_x
end
ũ_x = compute_u()[:]
pdf = kde(ũ_x)

plot(pdf.x, pdf.density)

gr()

using MeshCat
using CoordinateTransformations
using GeometryTypes: GeometryTypes, HyperRectangle, Vec, Point,
    HomogenousMesh, SignedDistanceField, HyperSphere, GLUVMesh
using Colors: RGBA, RGB
using MeshIO
using FileIO
QQ = [0.707107;-0.707107; -0.0; -0.0]
animate_traj(t_sim, samples)
