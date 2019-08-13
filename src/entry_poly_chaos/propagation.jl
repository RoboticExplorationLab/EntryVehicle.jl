#propagating dynamics for each sample

function propagation(fnc, x0, t_sim)
    t_sim, z_sim = ode78(fnc, x0, t_sim, points=:specified)
    n = 13 #number of states
    Z = zeros(n,length(t_sim))
    @show(t_sim)
    for i = 1:length(t_sim)
        Z[:,i] = z_sim[i]
    end
    return Z
end
