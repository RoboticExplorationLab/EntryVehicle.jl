# test unscented transform for dynamics propagation

using LinearAlgebra
using Plots
include("sig_point.jl")
include("unscented_prop.jl")
include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model_dis.jl")
include("traj_plots.jl")

n = 13 #number of states
Re = 3389.5
N = 500 #simulation length
dt = 0.01
T = 0.0:dt:N
store_state = zeros(length(T), n)
store_cov = zeros(n, n, length(T))


#initialization
x0 = [(3389.5+125)/Re; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 2.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal([1.0/Re; 1.0/Re; 1.0/Re; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001; 0.0001])

mu1 = x0
Q1 = Q0

function prop(mu1, Q1)
    store_state = zeros(length(T), n)
    store_cov = zeros(n, n, length(T))
    for i = 1:1:length(T)
        t = T[i]
        mu2, Q2 = Unscented_prop(mu1, Q1, t)
        store_state[i, :] = mu2'
        store_cov[:, :, i] = Q2
        mu1 = mu2
        Q1 = Q2
    end
    return store_state, store_cov
end

state, cov = prop(mu1, Q1)

plot(state[:, 1]*Re, state[:, 2]*Re)
