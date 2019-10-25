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
    Plots.plot(X[1, :]*Re, X[2, :]*Re, label="Spacecraft Trajectory")
    xlabel!("X [km]")
    ylabel!("Y [km]")
    title!("Spacecraft Trajectory in MCI - XY projection")
end

function plot_traj2(X)
    Plots.plot!(X[1, :]*Re, X[2, :]*Re, label="Spacecraft Trajectory", show=true)
    xlabel!("X [km]")
    ylabel!("Y [km]")
    title!("Spacecraft Trajectory in MCI - XY projection")
end

function plot_altitude(X, t_sim)
    Re = 3389.5
    R = zeros(length(t_sim))
    for i = 1:length(t_sim)
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-1)*Re
    end
    Plots.plot!(t_sim, R, label="Spacecraft Altitude", show=true)
    xlabel!("t [s]")
    ylabel!("altitude [km]")
    title!("Spacecraft Altitude")
end

function plot_quaternions(X)
    Plots.plot(X[4, :])
    Plots.plot!(X[5, :])
    Plots.plot!(X[6, :])
    Plots.plot!(X[7, :])
    xlabel!("t [s]")
    ylabel!("Quaternions")
    title!("Spacecraft Quaternions")
end

function plot_ang_vel(X)
    Plots.plot(X[11, :])
    Plots.plot!(X[12, :])
    Plots.plot!(X[13, :])
    xlabel!("t [s]")
    ylabel!("Angular Velocity Components")
    title!("Spacecraft Angular Velocity")
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
            V[i] = sqrt(X[8, i]^2+X[9, i]^2+X[10, i]^2)
        end
        Plots.plot(t_sim, V, label="Spacecraft Velocity")
        xlabel!("t [s]")
        ylabel!("Velocity [km/s]")
        title!("Spacecraft Velocity")
end

function vector_z_body(X, t_sim)
    Z_b = zeros(length(t_sim), 3)
    for i=1:length(t_sim)
        q = X[4:7, i]'
        Z_b[i, :] = qrot(qconj(q), [0.0;0.0;1.0])'
    end
    return(Z_b)
end
