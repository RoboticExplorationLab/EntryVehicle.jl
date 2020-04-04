using LinearAlgebra
using DifferentialEquations
using StaticArrays
using Plots
gr()

include("structures.jl")
include("dynamics\\entry_vehicle_3dof_energy.jl")
include("aerodynamics.jl")


geo_phoenix=SphereconeGeometry(70.0*pi/180, 1.325, 0.2)
phoenix=Vehicle("phoenix",geo_phoenix,600.0,SMatrix{3,3}([184.180 0.0 0.0;
                                        0.0 184.180 0.0;
                                        0.0 0.0 230.0]),
                    SMatrix{3,3}([0.00542947 0.0 0.0;
                                        0.0 0.00542947 0.0;
                                        0.0 0.0 0.00434783]),
                            SVector{3}([0.0,0.,0.42]))
table_CD, table_CL = drag_lift_table(geo_phoenix.δ, geo_phoenix.r_min,
                                    geo_phoenix.r_cone, phoenix.r_com)

α = 1:1:60
C = [table_CL[i]/table_CD[i] for i in α]
Plots.plot(α,C)

entry_params = (Vehicle=phoenix, A_ref=A_ref(geo_phoenix),L_ref=L_ref(geo_phoenix),
                    Env=Earth)

u = [0.0*pi/180]

g = 9.81
R0 = 6378.135*1e3

h0 = 121.9*1e3
r0 = h0+R0
V0 = 7623.5
e0 = e(r0/R0,V0/(sqrt(g*R0)))
hf = 7600
rf = hf+R0
Vf = 150.0
ef = e(rf/R0, Vf/sqrt(g*R0))
s0 = 1852*2035/R0 #radians

x0 = [r0/R0,0.0,0.0,-13*pi/180,0.0,s0]

espan = (e0, 0.989)
prob = ODEProblem(entry_vehicle_3dof_dynamics_energy!, x0, espan, entry_params)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

Plots.plot(sol, vars=5)


function σ(e,σ_0,σ,_f,e_0,e_f)
    #assume linear profile of σ with e for the guidance cycle
    return σ_0+((e-e_0)/(e_0-e_f))*(σ_f-σ_0)
end

function z(σ_0,s_obj)
    # compute the error at each step
    # s_obj is the downrange that we want to reach finally
    # σ is the current iterate
    espan = (e_0,e_f)
    prob = ODEProblem(entry_vehicle_3dof_dynamics_energy!,x0,espan,entry_params)
    sol = solve(prob)
    #EXTRACT s(f) = s(end)
    return (sf-s_obj)
end

function f(σ_0,s_obj)
    return 0.5*z(σ_0,s_obj)^2
end

function update_sigma(σ_0,σ_1,z_0,s_obj)
    λ = 1.0
    z_1 = z(σ_1,s_obj)
    f_1 = 0.5*z_1^2
    σ_2 = σ_1-λ*(z_1/(z_1-z_0))*(σ_1-σ_0)
    f_2 = f(σ_2)
    while f2>=f_1
        λ /= 2.0
        σ_2 = σ_1-λ*(z_1/(z_1-z_0))*(σ_1-σ_0)
        f_2 = f(σ_2)
    end
    return σ_2
end



################################################################################
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
        k1 = f(y_star, p)
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

# Test Traj NASA ###############################################################

using CSV
using DataFrames

file = "mars_edl_low_ld_v2.txt"
F = CSV.File(file)
df = DataFrame(F)
names(df)

Plots.plot(df.time, (sqrt.(df.xi.^2+df.yi.^2+df.zi.^2)))

Plots.plot(df.time, (sqrt.(df.vxi.^2+df.vyi.^2+df.vzi.^2)))

Plots.plot(df.time, df.vxi)
Plots.plot(df.time, df.vyi)
Plots.plot(df.time, df.vzi)
Plots.plot(df.time, df.mass) #constant
Plots.plot(df.time, df.piti)
Plots.plot(df.time, df.yawi)
Plots.plot(df.time, df.roli)

Plots.plot(df.time, df.sref) #constant

Plots.plot(df.time, df.pitr)
Plots.plot(df.time, df.rolr)
Plots.plot(df.time, df.yawr)

Plots.plot(df.time, df.alpha) #almost always constant here

Plots.plot(df.time, df.jdate) #constant makes sense

# 6 DOF simulation #############################################################

include("interpolation.jl")
include("dynamics\\quaternions.jl")
include("dynamics\\entry_vehicle_6dof.jl")

v_eci = [-1.6; 6.8; 0.0001]*1e3
β = acos((v_eci'*[0.0; 1.0; 0.0])/(norm(v_eci)))
x_b = [-cos(β);-sin(β); 0.0]
z_b = v_eci/(norm(v_eci))
y_b = [0.0; 0.0; 1.0]
M = hcat(x_b, y_b, z_b)
Q = mat2quat(M)
Q = qconj(Q) #quaternions initial

C_FX,C_FY,C_FZ,C_τX,C_τY,C_τZ,DX,DY,DZ = coeff_interp_6dof(geo_phoenix.δ,geo_phoenix.r_min,geo_phoenix.r_cone, phoenix.r_com)
coeff_interp=Aero_coeff_interp(C_FX,C_FY,C_FZ,C_τX,C_τY,C_τZ,DX,DY,DZ)

x0 = [(3389.5+125)*1e3; 0.0; 50.0; Q[1]; Q[2]; Q[3]; Q[4]; v_eci; 0.0; 0.0; 0.0]
tspan = (0.0,200.0)
entry_params_2 = (Vehicle=phoenix, A_ref=A_ref(geo_phoenix),L_ref=L_ref(geo_phoenix),
                    Env=Mars,coeff_interp=coeff_interp)
prob_2 = ODEProblem(entry_vehicle_6dof_dynamics!, x0, tspan, entry_params_2)
sol = solve(prob_2,Tsit5(),reltol=1e-8,abstol=1e-8)
Plots.plot(sol, vars=2)
