#File contains integration process

function integration(fnc, x0, Δt)
    Re = 3389.5
    #x0_bb = [6453.448873708138/Re; -0564.603617093604/Re; 0.0; 1.0; 0.0;0.0;0.0; 0.683661419343188; 7.814285780472696; 0.0; 0.1;0.0;0.0] #orbit 5 degrees before
    #x0_b = [6477.113353992618/Re; -0113.058434141366/Re; 0.0; 1.0;0.0;0.0;0.0;0.136899033611736;7.842940382904379; 0.0;0.0;0.0;0.0] #orbit 1 degree before so should land before
    t_sim, z_sim = ode78(fnc, x0, 0:0.1:Δt, points=:all) #or ODE45 otherwise. can change all or specified for the points we want.
    n = 13 #number of states
    Z = zeros(n,length(t_sim))

    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end


#Define Continuous Call Back when reaching the ground
function condition(u, t, integrator)
    norm(u[1:3])-1.0
end

function affect!(integrator)
    terminate!(integrator)
end

cb = ContinuousCallback(condition,affect!)

function integration2(fnc, x0, Δt)
    Re = 3389.5
    tspan = (0.0, Δt)
    prob = ODEProblem(fnc, x0, tspan)
    sol = DifferentialEquations.solve(prob, points=:specified, callback = cb) #save_everystep=false,dense=false)#reltol=1e-10,abstol=1e-10) #Tsit5(), Rodas4(), AutoTsit5(Rosenbrock23())
    n = 13 #number of states
    t_sim = sol.t
    z_sim = sol.u
    Z = zeros(n,length(t_sim))
    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return t_sim, Z
end
