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
N = 230 #simulation length
dt = 0.1
T = 0.0:dt:N
#T = vcat(0:0.1:230.0, 230.0:0.001:250.0, 250.1:0.0005:270)
store_state = zeros(length(T), n)
store_cov = zeros(n, n, length(T))
w = [0.0158*10^9; 0.0; 0.0; 0.0]
u = [0.0]

#initialization
θ = 91*pi/180 #rotation angle about z-axis body
M = [-sin(θ) cos(θ) 0.0;
     0.0 0.0 1.0;
     cos(θ) sin(θ) 0.0]
Q = mat2quat(M)
Q = qconj(Q)
x0 = [(3389.5+125)/Re; 0.0; 0.0; Q[1]; Q[2]; Q[3]; Q[4]; 0.0; 1.0; 0.0; 0.0; 0.0; 0.0]
Q0 = Diagonal([100.0/Re; 100.0/Re; 100.0/Re; 0.01; 0.01; 0.01; 0.01; 0.1; 0.01; 0.01; 0.001; 0.01; 0.01])

mu1 = x0
Q1 = Q0

function prop(mu1, Q1)
    store_state = zeros(n, length(T))
    store_cov = zeros(n, n, length(T))
    for i = 1:1:length(T)
        t = T[i]
        @show(t)
        if t<=230
            dt = 0.1
        elseif t <=250 && t>230
            dt = 0.001
        else
            dt = 0.0005
        end
        mu2, Q2 = transform_unscented(mu1, Q1, t, dt, u, w)
        store_state[:, i] = mu2'
        store_cov[:, :, i] = Q2
        mu1 = mu2
        Q1 = Q2
    end
    return store_state, store_cov
end

state, cov = prop(mu1, Q1)

#Evaluate solution
plot_traj(state)
plot_altitude(state, T)
norm_quaternions(state, T)
plot_quaternions(state)


function plot_all_ellipse()
    for i = 1:100:length(T)
        plot_ellipse(state[i, 1:2]*Re, cov[1:2, 1:2, i]*Re, 0.99)
    end
end

plot_all_ellipse()

plot_ellipse(state[1, 1:2]*Re, cov[1:2, 1:2, 1]*Re, 0.99)
