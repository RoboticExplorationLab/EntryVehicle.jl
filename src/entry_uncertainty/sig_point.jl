#sigma points extraction for UKF

function extract_sig_point(mu, sig)
    n = length(mu);
    sig_point = zeros(2*n+1, n)
    sig_point[1, :] = mu'
    kappa = 0.5;
    alpha = 0.5;
    lambda = (alpha^2)*(n+kappa)-n;
    C = cholesky(sig)
    mat = (C.U)
    for i=1:1:n
        sig_point[i+1, :] = mu' + sqrt(n+lambda)*mat[:, i]';
    end
    for i=n+1:1:2*n
        sig_point[i+1, :] = mu' - sqrt(n+lambda)*mat[:,i-n]';
    end
    return sig_point
end
