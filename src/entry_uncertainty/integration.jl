#File contains integration process

function integration(fnc, x0, Δt)
    Re = 3389.5
    #x0_bb = [6453.448873708138/Re; -0564.603617093604/Re; 0.0; 1.0; 0.0;0.0;0.0; 0.683661419343188; 7.814285780472696; 0.0; 0.1;0.0;0.0] #orbit 5 degrees before
    #x0_b = [6477.113353992618/Re; -0113.058434141366/Re; 0.0; 1.0;0.0;0.0;0.0;0.136899033611736;7.842940382904379; 0.0;0.0;0.0;0.0] #orbit 1 degree before so should land before
    t_sim, z_sim = ode78(fnc, x0, 0:0.1:Δt, points=:specified) #or ODE45 otherwise
    n = 13 #number of states
    Z = zeros(n,length(t_sim))

    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end

function integration2(fnc, x0, Δt)
    Re = 3389.5
    tspan = (0.0, Δt)
    prob = ODEProblem(fnc, x0, tspan)
    sol = DifferentialEquations.solve(prob) #reltol=1e-10,abstol=1e-10, verbose = false, points=:specified, saveat = 1.0) #Tsit5(), Rodas4(), AutoTsit5(Rosenbrock23())
    n = 13 #number of states
    t_sim = sol.t
    z_sim = sol.u
    Z = zeros(n,length(t_sim))
    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
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
