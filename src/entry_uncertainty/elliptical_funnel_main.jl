#Main file for propagation using UT and elliptical funnel
using Convex
using SCS
using Mosek
using Roots
using LinearAlgebra
using Plots
using Random
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
x0 = [0.0; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0])
#Diagonal([0.1/Re; 0.1/Re; 0.1/Re; 0.01; 0.01; 0.01; 0.01; 0.0001; 0.0001; 0.0001; 0.01; 0.01; 0.01])

#x0 = [2.0; 1.0]
#Q0 = Diagonal([1.0, 4.0])
Q0 = Matrix(Q0)

A1 = inv(sqrt(Q0))
b1 = -A1*x0

Δt = 100 #length simulation
dt = 1
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
        #X2 = prop_points(X1, t, dt, u, w)
        #X2 = X1
        A2, b2 = points2ellipse(X1)
        blist[:, i] = b1
        Alist[:, :, i] = A1
        centerlist[:, i] = -inv(A1)*b1
        A1 = A2
        b1 = b2
        XX[:, :, i] = X1
    end
    return Alist, blist, centerlist, XX
end

Alist, blist, centerlist, XX = propagation(A1, b1) #get conservative elliptical funnels for the trajectory

plot_traj_center(centerlist)

centerlist[:, 1]
centerlist[:, 100]

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

j = 40
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
    B[:, i] = [cos(angles[i]) - blist[1, j], sin(angles[i]) - blist[2, j]]
end

ellipse  = Alist[1:2, 1:2, j] \ B
plot!(ellipse[1, :], ellipse[2, :])
scatter!([centerlist[1, j]],[centerlist[2, j]] )
scatter!(XX[1, :, j], XX[2, :, j])

#test points extraction (ellipse2points)
A = Alist[1:2, 1:2, 30]
b = blist[1:2, 30]
M = -inv(A)*b
T = A'*A
W = eigvecs(T)
V = eigvals(T)
F = eigen(T)
W = F.vectors
V = F.values
v1 = M+(1/sqrt(V[1]))*W[:, 1]
v2 = M+(1/sqrt(V[2]))*W[:, 2]
v3 = M-(1/sqrt(V[1]))*W[:, 1]
v4 = M-(1/sqrt(V[2]))*W[:, 2]
scatter!([v1[1], v2[1], v3[1], v4[1]], [v1[2], v2[2], v3[2], v4[2]])


#test funnel in 2D only

Q0 = Diagonal([10.0/Re; 20.0/Re])
A1 = inv(sqrt(Q0))
b1 = -A1*x0[1:2]
n = length(x0[1:2])

Alist, blist, centerlist, XX = propagation(A1, b1)

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


function plot_traj_center(centerlist)
    X = zeros(length(T))
    Y = zeros(length(T))
    for j=1:length(T)
        X[j] = centerlist[1, j]#*Re
        Y[j] = centerlist[2, j]#*Re
    end
    plot(X, Y)
end

#=    R = W*Diagonal(sqrt.(z))*inv(W)
    @show(z)
    @show(b)
    for i = 1:n
        #going from dim 12 (ellipsoid) to dim 13 (propagation)
        points[1:3, 2*i-1] = C[1:3] + R[1:3, i] #position
        points[4:7, 2*i-1] = qmult(C[4:7], [sqrt(abs(1-R[4:6, i]'*R[4:6, i])); R[4:6, i]]) #quat (eigenvalues ?)
        points[8:13, 2*i-1] = C[8:13] + R[7:12, i] #remaining

        points[1:3, 2*i] = C[1:3] -R[1:3, i] #position
        points[4:7, 2*i] = qmult(C[4:7], [sqrt(abs(1-R[4:6, i]'*R[4:6, i])); -R[4:6, i]]) #quat (eigenvalues ?)
        points[8:13, 2*i] = C[8:13] - R[7:12, i] #remaining
    end=#
