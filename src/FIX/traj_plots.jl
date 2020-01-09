#File contains plot functions for entry trajectory
#2D Trajectory in Perifocal
#Altitude
#Entry Velocity
#Angular Velocity
#Quaternions
#Norm Quaternions
#Z body axis components in ECI

#Plot functions

function plot_traj(X)
    Plots.plot(X[1, :]*1e-3, X[2, :]*1e-3, legend = false)
    xlabel!("X [km]")
    ylabel!("Y [km]")
end

function plot_altitude(X, t_sim)
    R = zeros(length(t_sim))
    for i = 1:length(t_sim)
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-(3389.5*1e3))
    end
    Plots.plot(t_sim, R*1e-3, legend=false)
    xlabel!("Time [s]")
    ylabel!("Altitude [km]")
end

function plot_attack_angle(X, t_sim)
    n = length(t_sim)
    α = zeros(n)
    for i =1:1:n
        q = X[4:7, i]
        v = X[8:10, i]
        r = X[1:3, i]
        ω_mars = [0.0; 0.0; 7.095*10^(-5)]
        ω_mars = [0.0; 0.0; 0.0]
        v_rel = (v-cross(ω_mars, r)) #velocity of spacecraft wrt atm
        v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame
        α[i] = acos(v_body[3]/norm(v_body)) #radians
        α[i] = α[i]*180/pi
    end
    Plots.plot(t_sim, α)
end

function plot_total_attack_angle(X, t_sim)
    n = length(t_sim)
    α_t = zeros(n)
    for i =1:1:n
        q = X[4:7, i]
        v = X[8:10, i]
        r = X[1:3, i]
        ω_mars = [0.0; 0.0; 7.095*10^(-5)]
        #ω_mars = [0.0; 0.0; 0.0]
        v_rel = (v-cross(ω_mars, r)) #velocity of spacecraft wrt atm
        v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame
        α = acos(v_body[3]/norm(v_body))
        α_t[i] = α*180/pi
    end
    Plots.plot(t_sim, α_t, legend = false)
    ylabel!("Total AoA [°]")
    xlabel!("Time [s]")
end

function plot_entry_profile(X, t_sim)
    #defined as altitude vs relative velocity
    n = length(t_sim)
    R = zeros(length(t_sim))
    V = zeros(length(t_sim)) #relative vel
    ω_mars = [0.0; 0.0; 7.095*10^(-5)]
    #ω_mars = [0.0; 0.0; 0.0]
    for i = 1:1:n
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-(3389.5*1e3))
        q = X[4:7, i]
        v = X[8:10, i]
        r = X[1:3, i]
        v_rel = (v-cross(ω_mars, r)) #velocity of spacecraft wrt atm
        V[i] = norm(v_rel)
    end
    Plots.plot(V*1e-3, R*1e-3, legend=false)
    xlabel!("Relative velocity [km.s-1]")
    ylabel!("Altitude [km]")
end



function plot_quaternions(X, t_sim)
    Plots.plot(t_sim, X[4, :], label= "q0")
    Plots.plot!(t_sim, X[5, :], label= "q1")
    Plots.plot!(t_sim, X[6, :], label= "q2")
    Plots.plot!(t_sim, X[7, :], label= "q3")
    xlabel!("Time [s]")
    ylabel!("Quaternions")
end

function plot_ang_vel(X, t_sim)
    Plots.plot(t_sim, X[11, :], label = "omega_x")
    Plots.plot!(t_sim, X[12, :], label = "omega_y")
    Plots.plot!(t_sim, X[13, :], label = "omega_z")
    xlabel!("Time [s]")
    ylabel!("Angular Velocity Components")
end

function norm_quaternions(X, t_sim)
    N = zeros(length(t_sim))
    for i = 1:length(t_sim)
        N[i] = norm(X[4:7, i])
    end
    Plots.plot(t_sim, N)
end

function plot_vel(X, t_sim)
        V = zeros(length(t_sim))

        for i = 1:length(t_sim)
            V[i] = sqrt(X[8, i]^2+X[9, i]^2+X[10, i]^2)*1e-3
        end
        Plots.plot(t_sim, V, legend = false)
        xlabel!("Time [s]")
        ylabel!("Velocity [km/s]")
end

function vector_z_body(X, t_sim)
    Z_b = zeros(length(t_sim), 3)
    for i=1:length(t_sim)
        q = X[4:7, i]'
        Z_b[i, :] = qrot(qconj(q), [0.0;0.0;1.0])'
    end
    return(Z_b)
end

function plot_temp(X, t_sim)
    plot(t_sim, Z[14, :], legend = false)
    xlabel!("Time [s]")
    ylabel!("Temperature [K]")
end

function plot_acc(X, t_sim)
    Re = 3389.5
    R = zeros(length(t_sim))
    for i = 1:length(t_sim)
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-1)*Re
    end
    plot(Z[15, :], R, legend = false)
    ylabel!("Altitude [km]")
    xlabel!("Acceleration [g]")
end

function plot_mach_number(X, t_sim)
    Re = 3389.5*1e3
    n = length(t_sim)
    M = [norm(X[8:10, i])/(speed_sound(norm(X[1:3, i])-Re)) for i=1:n]
    Plots.plot(t_sim, M)
    xlabel!("Time [s]")
    ylabel!("Mach number")
end

function plot_mach_number_altitude(X, t_sim)
    Re = 3389.5*1e3
    n = length(t_sim)
    M = [norm(X[8:10, i])/(speed_sound(norm(X[1:3, i])-Re)) for i=1:n]
    A = [(norm(X[1:3, i])-Re)/(1e3) for i=1:n]
    Plots.plot(A, M, legend =false)
    xlabel!("Altitude [km]")
    ylabel!("Mach number")
end

function specific_energy(V, r)
    μ = 4.282837*1e13
    E = 0.5*V^2-μ/r
    return E
end

function plot_specific_energy(X, t_sim)
    V = [norm(X[8:10,i]) for i=1:1:length(t_sim)]
    R = [norm(X[1:3,i]) for i=1:1:length(t_sim)]
    EE = [specific_energy(V[i], R[i]) for i=1:1:length(t_sim)]
    Plots.plot(t_sim, EE)
    xlabel!("Time [s]")
    ylabel!("Specific energy [J]")
end
