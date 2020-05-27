# Integration Files containing all integrators

T = Float64
# Regular Fixed Step RK4 integrator
function rk4(f, y_0::Vector{T}, model::AbstractModel, dt::T, t_span::Vector{T})
    u = [0.0]   # TO BE MODIFIED LATER
    span = t_span[1]:dt:t_span[end]
    y = zeros(length(span), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(span)-1
        t = span[i]
        y_star = y[i, :]
        k1 = f(t, y_star, u, model)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, u, model)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, u, model)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, u, model)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return span, y'
end


function rk4(f, y_0, p, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    if length(y_0) == 1
        y[1, :] = [y_0]
    else
        y[1, :] = y_0
    end
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(t, y_star, p)   # Check that line for time
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(t+dt/2, y1, p)
        y2 = y_star+k2*dt/2
        k3 = f(t+dt/2, y2, p)
        y3 = y_star+k3*dt
        k4 = f(t+dt, y3, p)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y'
end


#=
# Unit Test on Duffing Oscillator
using LinearAlgebra
using Plots

# Simple Dynamics Function
function duffing(t,u,p)
    # p = α, β, δ, γ, ω
    du = zeros(2)
    α = -1.0 #p[1]
    β = 1.0 #p[2]
    δ = 0.2 #p[3]
    γ = 0.1 #p[4]
    ω = 1.0 #p[5]
    du[1] = u[2]
    du[2] =  γ*cos(t)-δ*u[2]-α*u[1]-β*u[1]^3 -0.01*(u[1]-1.0)-0.01*(u[2]-0.0) #+F(t)[1]
    return du
end

# Plot Trajetories with respect to time
function plot_trajectories_duffing!(T_sim, U_sim)
    plt = plot()
    Plots.plot!(T_sim, [U_sim[1, :], U_sim[2, :]],
                            xlabel = "Time",
                            ylabel = "Values",
                            label = ["position" "velocity"],
                            title = "Duffing Oscillator")
    display(plt)
    return nothing
end

# Plot Phase Portrait
function plot_phase_portrait_duffing!(U_sim)
    plt = plot()
    Plots.plot!(U_sim[1, :], U_sim[2, :],
                        xlabel = "position",
                        ylabel = "velocity",
                        label = "phase space trajectory, initial = $u_0",
                        title  = "Phase Portrait Duffing")
    display(plt)
    return nothing
end

# Actual Test
u_0 = [0.1, 10.0]
dt = 1e-2
t_span = [0.0, 50.0]
T_sim, U_sim = rk4(duffing, u_0, [0.0], dt, t_span)
plot_trajectories_duffing!(T_sim, U_sim)
plot_phase_portrait_duffing!(U_sim)

=#
