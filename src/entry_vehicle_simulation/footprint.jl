function footprint()
    state_end = zeros(36, 13)
    pos_end = zeros(36, 3)
    list = 0:10:350
    for i in 1:1:length(list)
        θ = list[i]*pi/180
        if i != 10
        #=M = [0.0 0.0 -1.0;
            -sin(θ) cos(θ) 0.0;
             cos(θ) sin(θ) 0.0] #Z body axis towards the planet =#
        M = [-sin(θ) cos(θ) 0.0;
             0.0 0.0 1;
             cos(θ) sin(θ) 0.0]  #Z body axis in the direction of the initial velocity
        Q = mat2quat(M)
        Q = qconj(Q)
        @show(Q)
        x0 = [(3389.5+125)/Re, 0.0/Re, 0.0, Q[1], Q[2], Q[3], Q[4], 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        t_sim, Z= integration(dyn, x0, Δt)
        state_end[i, :] = Z[:, end]'
        pos_end[i, :] = Z[1:3, end]'
        @show(i)
    end
    end
    return pos_end, state_end
end
