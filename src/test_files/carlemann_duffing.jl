#Carlemann linearization

function S(i, p)
    S = 0
    for k =0:1:i-1
        S += (p+1)-k
    end
    return S
end

function map_index(i, j, p)
    if i == 0
        return (j+1)
    else
        return S(i, p) + (j+1)
    end
end

a=1

function carlemann_matrix(p, t)
    #p is the order at which we truncate the linearization (polynomial order)
    n = Int64((p+1)*(p+2)/2)
    A = zeros(n, n)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p) #number of the line we are working at
            if k <= n && k >=1
                A[k, k] = -δ*j
                k1 = map_index(i+1, j-1, p)
                if k1 <= n && k1 >=1
                    A[k, k1] = -j*α
                end
                k2 = map_index(i+3, j-1, p)
                if k2 <= n && k2 >=1
                    A[k, k2] = -j*β
                end
                k3 = map_index(i, j-1, p)
                if k3 <= n && k3 >=1
                    A[k, k3] = -j*γ*cos(t)
                end
                k4 = map_index(i-1, j+1, p)
                if k4 <= n && k4 >=1
                    A[k, k4] = i
                end
            end
        end
    end
    return A
end

α = -1.0
β = 1.0
γ = 0.2
δ = 0.1
A = carlemann_matrix(2, pi)

function basis(x, p)
    n = Int64((p+1)*(p+2)/2)
    base = zeros(n)
    for i = 0:1:p
        for j = 0:1:p-i
            k = map_index(i, j, p)
            base[k] = (x[1]^i)*(x[2]^j)
        end
    end
    return base
end

function integration_carlemann(x_ini, t_ini, t_end, dt)
    T = t_ini:dt:t_end
    X = zeros(length(x_ini), length(T))
    X[:, 1] = x_ini
    for q=1:1:length(T)-1
        A = carlemann_matrix(p, T[q])
        Phi = exp(A*dt)
        k1 = map_index(1, 0, p)
        k2 = map_index(0, 1, p)
        X[:, q+1] = hcat(Phi[k1, :], Phi[k2, :])'*basis(X[:, q], p)
    end
    return X
end

t_ini = 0.0
t_end = 100.0
dt = 0.01
p = 6
x_ini = [0.1 10.0]
X = integration_carlemann(x_ini, t_ini, t_end, dt)

Plots.plot(X[1, :], X[2, :])
