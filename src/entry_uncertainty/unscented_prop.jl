#test Unscented Propagation on entry model

function Unscented_prop(mu_1, sig_1, t, dt, w)
    #Parameters
    n_Q, m_Q = size(sig_1)
    #compute sigma points
    n = length(mu_1);
    sig_point_0 = zeros(2*n+1, n)
    sig_point_0[1, :] = mu_1'
    alpha = 0.5;
    kappa = 17;
    beta = 5;
    lambda = (alpha^2)*(n+kappa)-n;
    weight_mean = zeros(2*n+1)
    weight_cov = zeros(2*n+1)
    weight_mean[1] = lambda/(n+lambda);
    weight_cov[1] = (lambda/(n+lambda))+ (1-alpha^2+beta)
    C = cholesky(sig_1)
    mat = C.U
    for i=1:1:n
        sig_point_0[i+1, :] = mu_1' + sqrt(n+lambda)*mat[:, i]';
        weight_mean[i+1] = 1/(2*(n+lambda));
        weight_cov[i+1] = 1/(2*(n+lambda));
    end
    for i=n+1:1:2*n
        sig_point_0[i+1, :] = mu_1' - sqrt(n+lambda)*mat[:,i-n]';
        weight_mean[i+1] = 1/(2*(n+lambda));
        weight_cov[i+1] = 1/(2*(n+lambda));
    end
    sig_point_1 = zeros(2*n+1, n_Q);
    for i=1:1:2*n+1
        sig_point_1[i, :] = dyna_dis(t, dt, sig_point_0[i,:]', [0.0], w);
    end
    mu_bar = zeros(n_Q)
    sig_bar = zeros(n_Q, n_Q);
    for i=0:1:2*n
        mu_bar = mu_bar + weight_mean[i+1]*sig_point_1[i+1, :];
    end
    for i=0:1:2*n
        sig_bar = sig_bar + weight_cov[i+1]*((sig_point_1[i+1, :]-mu_bar)*(sig_point_1[i+1, :]-mu_bar)');
    end
    return mu_bar, sig_bar
end

function Unscented_prop_simple(mu_1, sig_1, t, w)
    n = length(mu_1)
    α = 0.01
    mu_2 = dyna_dis(t, dt, mu_1, [0.0], w)
    C = cholesky(sig_1)
    mat = C.U
    M = zeros(size(mat))
    for i = 1:1:n
        M[:, i] = dyna_dis(t, dt, mu_1+α*mat[:, i], [0.0], w)
        M[:, i] = (M[:, i] - mu_2)/α
    end
    sig_2 = M*M'
    return mu_2, sig_2
end

function Unscented_prop_simple2(mu_1, sig_1, t, w)
    n = length(mu_1)
    sig_point_1 = extract_sig_point(mu_1, sig_1)
    #lines of sig_point_1 contains the sigma points
    sig_point_2 = zeros(size(sig_point_1))
    for i = 1:1:2*n+1
        sig_point_2[i, :] = dyna_dis(t, dt, sig_point_1[i, :]', [0.0], w)
    end
    kappa = 17
    alpha = zeros(2*n+1)
    alpha[1] = kappa/(kappa+n)
    for i = 1:1:2*n
        alpha[i+1] = 0.5*(1/(n+kappa))
    end
    mu_2 = zeros(n)
    sig_2 = zeros(size(sig_1))
    for i = 1:1:2*n+1
        mu_2 = mu_2 + alpha[i]*sig_point_2[i, :]
    end
    for i = 1:1:2*n+1
        sig_2 = sig_2 + alpha[i]*(sig_point_2[i, :]-mu_2)*(sig_point_2[i, :]-mu_2)'
    end
    return mu_2, sig_2
end

function transform_unscented(mu_1, sig_1, t, dt, u, w)
    n = length(mu_1)
    α = 2.5
    #construct sigma_points
    sigma_points = zeros(n, 2*n+1) #columns are sigma points
    sigma_points[:, 1] = mu_1
    L = cholesky(sig_1).L
    for i = 1:1:n
        sigma_points[:, i+1] = mu_1 + α*L[:, i]
    end
    for i = n+1:1:2*n
        sigma_points[:, i+1] = mu_1 - α*L[:, i-n]
    end
    #get the mean
    mu_2 = zeros(n)
    for i=1:1:2*n+1
        mu_2 += dyna_dis(t, dt, sigma_points[:, i], u, w)/(2*n+1)
    end
    #get the covariance
    sig_2 = zeros(size(sig_1))
    for i=1:1:2*n+1
        sig_2 += (dyna_dis(t, dt, sigma_points[:, i], u, w)-mu_2)*(dyna_dis(t, dt, sigma_points[:, i], u, w)-mu_2)'/(2*n+1)
    end
    return mu_2, sig_2
end
