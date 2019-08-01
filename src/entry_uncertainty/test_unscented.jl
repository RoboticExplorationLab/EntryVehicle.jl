# test unscented transform for dynamics propagation

using LinearAlgebra
using Plots

include("sig_point.jl")
include("unscented_prop.jl")
include("quaternions.jl")
include("aerodynamic_coeff.jl")
include("entry_model_dis.jl")
include("traj_plots.jl")
include("plot_ellipse.jl")

n = 13 #number of states
Re = 3389.5
N = 450 #simulation length
dt = 0.1
T = vcat(0.0:0.1:332, 333:0.0001:450)
store_state = zeros(length(T), n)
store_cov = zeros(n, n, length(T))
w = [0.0158*10^9; 0.0; 0.0; 0.0]

#initialization
x0 = [(3389.5+125)/Re; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0; 2.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal([100.0/Re; 100.0/Re; 100.0/Re; 0.0001; 0.0001; 0.0001; 0.0001; 0.1; 0.1; 0.1; 0.0001; 0.0001; 0.0001])

mu1 = x0
Q1 = Q0

function prop(mu1, Q1)
    store_state = zeros(length(T), n)
    store_cov = zeros(n, n, length(T))
    for i = 1:1:length(T)
        t = T[i]
        if t<=332
            dt = 0.1
        elseif t>332
            dt = 0.0001
        end
        @show(t)
        mu2, Q2 = Unscented_prop(mu1, Q1, t, dt, w)
        store_state[i, :] = mu2'
        store_cov[:, :, i] = Q2
        mu1 = mu2
        Q1 = Q2
    end
    return store_state, store_cov
end

state, cov = prop(mu1, Q1)

Plots.plot(state[:, 1]*Re, state[:, 2]*Re)

function plot_all_ellipse()
    for i = 1:100:length(T)
        plot_ellipse(state[i, 1:2]*Re, cov[1:2, 1:2, i]*Re, 0.99)
    end
end

plot_all_ellipse()

plot_ellipse(state[1, 1:2]*Re, cov[1:2, 1:2, 1]*Re, 0.99)
