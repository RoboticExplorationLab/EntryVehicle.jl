using LinearAlgebra
using DifferentialEquations
using Plots
gr()

include("quaternions.jl")


function sys(u, p, t)
    du = zeros(7) 
    I1 = p[1] #assume principal axes
    I2 = p[2]
    I3 = p[3]
    I4 = p[4]
    I = Diagnoal([I1; I2; I3])
    I_inv = Diagnoal([1/I1; 1/I2; 1/I3])
    ω = u[5:7] #angular velocity at current step
    q = u[1:4] #quaternion at current step
    du[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7]))
    du[1:4] = 0.5*qmult(q, [0; ω])
    return du
end


function rk4(f, y_0, dt, t_span)
    T = t_span[1]:dt:t_span[end]
    y = zeros(length(T), length(y_0))
    y[1, :] = [y_0]
    for i=1:1:length(T)-1
        t = T[i]
        y_star = y[i, :]
        k1 = f(y_star, t)
        y1 = y_star+k1*dt/2 #intermediate evaluation value
        k2 = f(y1, t+dt/2)
        y2 = y_star+k2*dt/2
        k3 = f(y2, t+dt/2)
        y3 = y_star+k3*dt
        k4 = f(y3, t+dt)
        m = (k1+2*k2+2*k3+k4)/6 #slope average
        y[i+1, :] = y_star + m*dt
    end
    return T, y
end

function f(y, t)
    return -2*y
end

t_span = [0.0; 5.0]
y_0 = 3.0
dt = 0.6
t_sim, y = rk4(f, y_0, dt, t_span)
plot(t_sim, y)
T = 0.0:0.01:5.0
plot!(T, [3.0*exp(-2*t) for t in T])=#
