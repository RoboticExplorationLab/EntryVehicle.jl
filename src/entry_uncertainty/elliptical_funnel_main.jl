#Main file for propagation using UT and elliptical funnel
using Convex
using SCS
using Roots
using LinearAlgebra

include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model_dis.jl") #use discrete models
include("ellipse_points.jl") #used for extracting points and building ellipses
include("prop_points.jl")

#Initialization
Re = 3389.5
θ = 91*pi/180 #rotation angle about z-axis
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re, 0.0, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
Q0 = Diagonal([10.0/Re; 10.0/Re; 10.0/Re; 0.01; 0.01; 0.01; 0.01; 0.1; 0.01; 0.01; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = Q0
b1 = -inv(A1)*x0

Δt = 0.5 #length simulation
dt = 0.1
T = 0.0:dt:Δt

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

n = length(x0)

function propagation(A1, b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1) #return a set of points
        @show(size(X1))
        X2 = prop_points(X1, t, dt, u, w)
        A2, b2 = points2ellipse(X2)
        @show(size(b2))
        blist[:, i] = b2
        Alist[:, :, i] = A2
        A1 = A2
        b1 = b2
    end
    return Alist, blist
end

Alist, blist = propagation(A1, b1) #get conservative elliptical funnels for the trajectory
