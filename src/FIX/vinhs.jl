#Vinh's model for entry vehicle
#simplified version considering:
#-The oblateness of the planet is very small
#-Transport terms due to rotation between frames are small
#STANDARD SIMPLIFIED EQUATIONS OF MOTION
#See Advances in Spacecraft Atmospheric Entry Guidance (paper)

using LinearAlgebra
using Plots

include("aero.jl")
include("quaternions.jl")

a = 1

function vinhs_model(t, u, p)
    du = zeros(length(u))
    r, θ, ϕ, v, γ, ψ = u #state vector
    #theta is basically longitude
    #phi is the latitude
    σ = p[1] #bank angle (control input)


    @show(t)

    ω = 7.095*10^(-5) #rad.s-1 MARS
    μ =  4.282837*1e13 #m3.s-2

    g = μ/(r^2)
    C_γ = 2*ω*cos(ψ)*cos(ϕ) #Coriolis
    C_ψ = 2*ω*(tan(γ)*sin(ψ)*cos(ϕ)-sin(ϕ)) #Coriolis

    C_D = table_CD[16]
    C_L = table_CL[16]

    m = 600.0 # kg
    Re = 3389.5*1e3 #m
    h = r - Re #m
    D = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_D
    L = 0.5*exponential_atmosphere(h)*v^2*(A_ref/m)*C_L

    du[1] = v*sin(γ)
    du[2] = (v/r)*(cos(γ)*cos(ψ)/cos(ϕ))
    du[3] = (v/r)*cos(γ)*sin(ψ)
    du[4] = -D-g*sin(γ)
    du[5] = (1/v)*(L*cos(σ)-(g-(v^2)/r)*cos(γ)) + C_γ
    du[6] = -(1/(v*cos(γ)))*(L*sin(σ)+(v^2/r)*(cos(γ)^2)*cos(ψ)*tan(ϕ))+C_ψ
    if h>0.0
        return du
    else
        return zeros(6)
    end
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
        k1 = f(t, y_star, p)
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


δ = 70.0*pi/180.0
r_min = 0.2
r_cone = 0.3
r_G = [0.2; 0.0; 0.3]
A_ref = pi*r_cone^2
table_CD, table_CL = drag_lift_table(δ, r_min, r_cone, r_G)

u0 = [(122+3389.5)*1e3; 0.0; 0.0; 11.032*1e3; -15.0*pi/180; 0.0]
t_sim4, Z4 = rk4(vinhs_model,u0, [45.0*pi/180], 0.01, [0.0; 200.0])

Plots.plot(Z4[1, :])
Plots.plot(Z4[2, :])
Plots.plot(Z4[3, :])
Plots.plot(Z4[4, :])
Plots.plot(Z4[5, :])
Plots.plot(Z4[6, :])
