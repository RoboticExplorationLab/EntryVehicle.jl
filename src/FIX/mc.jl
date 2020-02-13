using Random
using Distributions

###############################################################################
######################### Monte Carlo Simulation ##############################
###############################################################################

#Gaussian here only so far for the paper
function generate_samples(x_0, Q, M)
    #M number of samples
    n = length(x_0)
    #rng = MersenneTwister(1234)
    #MVN = MvNormal(x_0, Q)
    X_samples = zeros(n, M)
    #univariate_D_vector = [Uniform(x_0[i]-sqrt(Q[i,i]),x_0[i]+sqrt(Q[i,i])) for i=1:length(x_0)]
    #D = Product(univariate_D_vector)
    MVN = MvNormal(x_0, Q/9)
    #=c = 0
    while c < M
        @show(c)
        X = rand!(D, zeros(n))
        if (X-x_0)'*inv(Q)*(X-x_0) <= 1.0
            X_samples[:, c+1] = X
            c +=1
        end
    end =#
    X_samples = rand!(MVN, X_samples)
    return X_samples
end

function prop_MC(X_samples, t_start, t_end, dt, p, model)
    n, M = size(X_samples)
    #saveAT = 1.0
    T_mc = t_start:dt:t_end
    traj = zeros(n, length(T_mc), M)
    for i=1:1:M
        @show(i)
        x_ini = X_samples[:, i]
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(model, x_ini, p, dt, [t_start, t_end])
        traj[:, :, i] = Z
    end
    return traj, T_mc
end

function mean_var_MC(traj)
    n, t, M = size(traj)
    avg = zeros(n, t)
    var = zeros(n, n, t)
    for i=1:1:t
        @show(i)
        S = zeros(n)
        V = zeros(n, n)
        for j=1:1:M
            S += traj[:, i, j]
        end
        avg[:, i] = S/M
        for k =1:1:M
            V+= (traj[:, i, k]-avg[:, i])*(traj[:, i, k]-avg[:, i])'
        end
        var[:, :, i] = V/M
    end
    return avg, var
end


function simu_MC(x0, Q0, M, model, p, t0, tf, dt_mc)
    X_samples = generate_samples(x0, Q0, M)
    traj, T_mc = prop_MC(X_samples, t0, tf, dt_mc, p, model)#counting prop
    avg, var, V = mean_var_MC(traj)
    return T_mc, traj, avg, var
end

function sig(var)
    #return 3*sig with respect to time for each state variable
    n, n, t = size(var)
    V = zeros(n, t)
    for i=1:1:n
        vec = [sqrt(var[i, i, j]) for j=1:1:t]
        V[i, :] = vec
    end
    return V
end


#= Example Vinhs_model
x0 = [(125+3389.5)*1e3; 0.0; 0.0; 7.032*1e3; -15.0*pi/180; 0.0]
Q0 = Diagonal([(50.0^2); (0.0017)^2; (0.0017)^2; (1.0)^2; (0.0001)^2; (0.0001)^2])
M = 100
model = vinhs_full_model
p = [45*pi/180]
t0 = 0.0
tf = 80.0
dt_mc = 0.01

T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M, model, p, t0, tf, dt_mc)

Plots.plot(T_MC, [traj_MC[1, :, i] for i=1:1:M])

V = [sqrt(var_MC[1, 1, i]) for i=1:1:length(T_MC)]
Plots.plot(T_MC, 3*V)
Plots.plot!(T_MC, -3*V) =#

#= Example Duffing oscillator

x0 = [0.1; 0.1]
Q0 = Diagonal([0.01^2; 0.01^2])
M = 1000
model = duffing
p = [-1.0; 1.0; 0.2; 0.1; 1.0]
t0 = 0.0
tf = 25.0
dt_mc = 0.01

T_MC, traj_MC, m_MC, var_MC = simu_MC(x0, Q0, M, model, p, t0, tf, dt_mc) =#


########################Additional functions 6 DOF##############################

function point12_13_mc(X)
    X13 = zeros(length(X)+1)
    X13[1:3] = X[1:3]
    X13[4:7] = euler2quat(X[4:6])
    X13[8:13] = X[7:12]
    return X13
end

function point13_12_mc(X)
    X12 = zeros(length(X)-1)
    X12[1:3] = X[1:3]
    X12[4:6] = quat2euler(X[4:7])
    X12[7:12] = X[8:13]
    return X12
end

function point13_12_mc_full(Z)
    n, m = size(Z)
    ZZ = zeros(n-1, m)
    for i=1:m
        ZZ[:, i] = point13_12_mc(Z[:, i])
    end
    return ZZ
end

function prop_MC_entry(X_samples, t_start, t_end, dt, p, model)
    n, M = size(X_samples)
    #saveAT = 1.0
    TT = t_start:dt:t_end
    traj = zeros(n, length(TT), M)
    for i=1:1:M
        @show(i)
        x_ini = point12_13_mc(X_samples[:, i])
        #prob = ODEProblem(duffing!,u0,tspan,M)
        #sol = DifferentialEquations.solve(prob, saveat = saveAT, abstol = 1e-9, reltol = 1e-9)
        t_sim, Z = rk4(model, x_ini, p, dt, [t_start, t_end])
        traj[:, :, i] = point13_12_mc_full(Z)
    end

    return traj, TT
end

#=
function mean_var_MC(traj)
    n, t, M = size(traj)
    avg = zeros(n, t)
    var = zeros(n, n, t)
    for i=1:1:t
        @show(i)
        S = zeros(n)
        V = zeros(n, n)
        for j=1:1:M
            S += traj[:, i, j]
        end
        avg[:, i] = S/M
        for k =1:1:M
            V+= (traj[:, i, k]-avg[:, i])*(traj[:, i, k]-avg[:, i])'
        end
        var[:, :, i] = V/M
    end
    return avg, var
end =#

function simu_mc_6dof(x0, Q0, M_mc, model, p, t0, tf, dt_mc)
    X_samples = generate_samples(x0, Q0, M_mc)
    traj_mc, T_mc = prop_MC_entry(X_samples, t0, tf, dt_mc, p, model)
    m_mc, var_mc = mean_var_MC(traj_mc)
    return T_mc, traj_mc, m_mc, var_mc
end
