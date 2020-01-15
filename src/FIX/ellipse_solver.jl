#home-made algo

function matrix(M, a, z)
    return (M*a-z)*a'+a*(M*a-z)'
end

function Hessian(M, z, A, t, u)
    n = length(z)
    m = length(u)
    grad_MM = 2*sum(u[i]*A[:, i]*A[:, i]' for i=1:1:m) - inv(M*M) #nxn
    grad_zM = -sum(u[i]*(A[:, i]+A[:, i]') for i=1:1:m) #
    grad_tM = zeros()
    grad_uM = [matrix(M, A[:, i], z) for i = 1:1:m] #n*n en ligne m fois
    grad_zz = sum(u)*Diagonaal(ones(n))
    grad_tz = zeros(n, m)
    grad_
    
