#extract points from ellipse
#need at least 2*n+1 points for an ellipsoid in n dimension ?

using Convex
using SCS

function ellipse2points(A, b)
    #A is the matrix we obtain from the previous step
    n = length(A[1, :])
    points = zeros(n, 2*n)
    M = -inv(A)*b #center of the ellipse considered
    P = inv(A'*A)
    W = eigvecs(P) #A symmetric positive definite
    z = eigvals(P)
    for i = 1:n
        #f(λ) = norm(λ*A*W[:, i], 2)-1
        #z = find_zeros(f, -10^9, 10^9) #refine the boundaries on this problem
        #@show(z)
        points[:, 2*i-1] = M + sqrt(z[i])*W[:, i]
        points[:, 2*i] = M - sqrt(z[i])*W[:, i]
    end
    return points #return 2n points so far
end

function points2ellipse(X)
    #X is the set of points propagated through the non linear dynamics
    n, m = size(X);
    #Define Convex Optimization Problem using Convex.jl
    A = Semidefinite(n)
    b = Variable(n)
    problem = maximize(logdet(A), vcat([norm(A*X[:, i]+b, 2)<=1 for i = 1:1:m], [A[1, 2]==A[2, 1]]))
    Convex.solve!(problem, SCSSolver())
    b = b.value
    A = A.value
    return A, b
end

#function test
#=
x = [1.55  2.0 1.0;
      2.25  2.35 1.0;
     2.2   2.2 1.0;
     2.25 2.1 1.0;
     2.0  2.3 1.0;
      2.4  2.2 1.0 ]';
A, b = points2ellipse(x)
points = ellipse2points(A, b)
A = A.value[1:2, 1:2]
b = b.value[1:2]
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end
ellipse  = A \ B
scatter(x[1,:], x[2,:])
plot!(ellipse[1, :], ellipse[2, :])
scatter!(points[1, :], points[2, :]) =#
