# Plot covariance ellipse

function plot_ellipse(mu, sigma, P)
    r = 1.0
    #r = sqrt(-2*log(1-P))
    alpha = 0:0.01:2*pi

    M = (cholesky(sigma).U)
    cxy = zeros(length(alpha), 2)
    for i = 1:1:length(alpha)
        cxy[i, :] = M*[r*cos(alpha[i]); r*sin(alpha[i])] + mu
    end
    Plots.plot(cxy[:, 1], cxy[:, 2])
end

plot_ellipse([0.0;0.0], [0.7 0.73; 0.73 1.1], 0.9999999999999999)
plot_ellipse([0.0; 0.0], 0.1*[0.7 0.73; 0.73 1.1], 0.96)

V = eigvecs([0.7 0.73; 0.73 1.1])
val = eigvals([0.7 0.73; 0.73 1.1])

V*Diagonal(val)*inv(V)
R = V*Diagonal(sqrt.(val))*inv(V)
R*R

scatter!([R[1, 1], R[2, 1]], [R[1, 2], R[2, 2]])
