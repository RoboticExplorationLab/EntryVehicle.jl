#Main file for propagation using UT and elliptical funnel
using Convex
using SCS
using Roots
using LinearAlgebra
using Plots
pyplot()

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
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal([10.0/Re; 20.0/Re; 30.0/Re; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01; 0.01])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0

Δt = 0.05 #length simulation
dt = 0.01
T = 0.0:dt:Δt

u = [0.0]
w = [0.0158*10^9; 0.0; 0.0; 0.0]

n = length(x0)

function propagation(A1, b1)
    blist = zeros(n, length(T))
    Alist = zeros(n, n, length(T))
    centerlist = zeros(n, length(T))
    XX = zeros(n, 2*n+1, length(T))
    for i=1:1:length(T)
        t = T[i]
        @show(t)
        X1 = ellipse2points(A1, b1) #return a set of points
        X2 = prop_points(X1, t, dt, u, w)
        A2, b2 = points2ellipse(X2)
        blist[:, i] = b2
        Alist[:, :, i] = A2
        centerlist[:, i] = -inv(A2)*b2
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX
end

Alist, blist, centerlist, XX = propagation(A1, b1) #get conservative elliptical funnels for the trajectory

centerlist[:, 1]
centerlist[:, end]

#test plot ellipse in 2D (positions)
function plott()
    for j=1:1:length(blist[1, :])
        angles = 0.0:0.01:2*pi
        B = zeros(2, length(angles))
        for i = 1:1:length(angles)
            B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
        end

        ellipse  = Alist[1:2, 1:2, j] \ B
        plot!(ellipse[1, :], ellipse[2, :])
    end
    return 0
end

j = 2
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
plot(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!([XX[1, :, j]], [XX[2, :, j]])
