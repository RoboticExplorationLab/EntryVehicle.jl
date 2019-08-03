# Plot covariance ellipse
using Plots
using LinearAlgebra
pyplot()

function plot_ellipse(mu, sigma, P)
    r = 1.0
    #r = sqrt(-2*log(1-P))
    alpha = 0:0.01:2*pi

    M = (cholesky(sigma).U)
    cxy = zeros(length(alpha), 2)
    for i = 1:1:length(alpha)
        cxy[i, :] = M*[r*cos(alpha[i]); r*sin(alpha[i])] + mu
    end
    Plots.plot!(cxy[:, 1], cxy[:, 2])
end

A = [2 5.0; 18.0 4]
b = [1.0; 1.0]
Q0 = Diagonal([100.0/Re; 100.0/Re; 100.0/Re; 0.01; 0.01; 0.01; 0.01; 0.1; 0.01; 0.01; 0.001; 0.01; 0.01])

angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end
ellipse  = A \ B
plot(ellipse[1, :], ellipse[2, :])

scatter!(-(inv(A)*b)[1, :], -(inv(A)*b)[1, :])

V = eigvecs(A)
val = eigvals(A)
V*Diagonal(val)*inv(V)
R = V*Diagonal(sqrt.(val))*inv(V)
scatter!([R[1, 1], R[2, 1]], [R[1, 2], R[2, 2]])

scatter!(V[1, :], V[2, :])

plot_ellipse([0.0;0.0], [0.7 0.73; 0.73 1.1], 0.9999999999999999)
plot_ellipse([0.0; 0.0], 0.1*[0.7 0.73; 0.73 1.1], 0.96)

V = eigvecs([0.7 0.73; 0.73 1.1])
val = eigvals([0.7 0.73; 0.73 1.1])

V*Diagonal(val)*inv(V)
R = V*Diagonal(sqrt.(val))*inv(V)
R*R

scatter!([R[1, 1], R[2, 1]], [R[1, 2], R[2, 2]])

AA = A[1, 1]
BB = 2*A[1, 2]
CC = A[2, 2]

aa = -sqrt(2*(4*AA*CC-BB^2)*((AA+CC)+sqrt((AA-CC)^2+BB^2)))/(BB^2-4*AA*CC)
bb = -sqrt(2*(4*AA*CC-BB^2)*((AA+CC)-sqrt((AA-CC)^2+BB^2)))/(BB^2-4*AA*CC)

if -det(A) < 0
    eta = 1
else
    eta = -1
end
e = sqrt((2*sqrt((AA-CC)^2+BB^2))/(eta*(AA+CC)+sqrt((AA-CC)^2+BB^2)))

scatter!(aa*(1-e)*V[1, :], aa*(1-e)*V[2, :])

using Roots

f(λ) = norm(λ*A*V[:, 2]+b)-1
find_zeros(f, -100, 100)

scatter!(V[1, :], V[2, :])

#test new form ellipse

A = [0.7 0.73/2 0.0; 0.73/2 1.1 0.0; 0.0 0.0 1]
cholesky(A)

V = eigvecs(A)
val = eigvals(A)
V*Diagonal(val)*inv(V)
R = V*Diagonal(sqrt.(val))*inv(V)

W = eigvecs(R)

#test covariance like matrix

D = [1.0 0.0;0.0 3.0]
b = [0.0; 0.0]
angles = 0.0:0.01:2*pi
B = zeros(2, length(angles))
for i = 1:1:length(angles)
      B[:, i] = [cos(angles[i]) - b[1], sin(angles[i]) - b[2]]
end
ellipse = inv(D) \ B
plot(ellipse[1, :], ellipse[2, :])

V = eigvecs(D)
val = eigvals(D)
V*Diagonal(val)*inv(V)
R = V*Diagonal(sqrt.(val))*inv(V)

W = eigvecs(R)
F = eigvecs(D)

scatter!(W[1, :], W[2, :])
scatter!(F[1, :], F[2, :])

using Roots

f(λ) = norm(λ*D*F[:, 2]+b)-1
z = find_zeros(f, -100, 100)

scatter!(z[1]*F[1, :], z[1]*F[2, :])
