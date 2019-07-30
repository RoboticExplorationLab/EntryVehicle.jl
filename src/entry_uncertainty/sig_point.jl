#sigma points extraction for UKF

function extract_sig_point(mu, sig)
    n = length(mu);
    sig_point = [mu'];
    alpha = 1;
    kappa = 1;
    lambda = (alpha^2)*(n+kappa)-n;
    C = cholesky(sig_1)
    mat = C.U
    for i=1:1:n
        sig_point[i+1, :] = mu' + sqrt(n+lambda)*mat[:, i]';
    end
    for i=n+1:1:2*n
        sig_point[i+1, :] = mu' - sqrt(n+lambda)*mat[:,i-n]';
    end
    return sig_point
end
