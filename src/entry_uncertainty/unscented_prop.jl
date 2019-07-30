#test Unscented Propagation on entry model

function Unscented_prop(mu_1, sig_1, t)
    #Parameters
    n_Q, m_Q = size(sig_1)
    #compute sigma points
    n = length(mu_1);
    sig_point_0 = zeros(2*n+1, n)
    sig_point_0[1, :] = mu_1'
    alpha = 0.7;
    kappa = 10;
    beta = 2;
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
        sig_point_1[i, :] = dyna_dis(t, dt, sig_point_0[i,:]', [0.0]);
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
