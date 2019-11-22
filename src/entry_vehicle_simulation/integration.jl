#File contains integration procedures

####################################
######Continuous Callback Code######
####################################
#Used to stop the simulation when
#reahcing ground

function condition(u, t, integrator)
    norm(u[1:3])-1.0
end

function affect!(integrator)
    terminate!(integrator)
end

cb = ContinuousCallback(condition,affect!)

####################################
######ODE.jl Integrator ############
####################################


function integration(fnc, x0, Δt)
    Re = 3389.5
    t_sim, z_sim = ode78(fnc, x0, 0:0.1:Δt, points=:all) #or ODE45 otherwise. can change all or specified for the points we want.
    n = 13 #number of states
    Z = zeros(n,length(t_sim))
    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end


####################################
###### DiffEq.jl Integrator ########
####################################


function integration2(fnc, x0, Δt)
    Re = 3389.5
    tspan = (0.0, Δt)
    prob = ODEProblem(fnc, x0, tspan)
    sol = DifferentialEquations.solve(prob, BS3(), points=:specified, callback = cb) #save_everystep=false,dense=false)#reltol=1e-10,abstol=1e-10) #Tsit5(), Rodas4(), AutoTsit5(Rosenbrock23())
    n = 13 #number of states
    t_sim = sol.t
    z_sim = sol.u
    Z = zeros(n,length(t_sim))
    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end


####################################
######Runge-Kutta 4 Integrator######
####################################


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
