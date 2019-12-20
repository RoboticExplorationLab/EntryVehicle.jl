#File contains plot functions for entry trajectory
#2D Trajectory in Perifocal
#Altitude
#Entry Velocity
#Angular Velocity
#Quaternions
#Norm Quaternions
#Z body axis components in ECI

#Plot functions

function plot_sphere2(X, t_sim)
    Re = 3389.5
    R = zeros(length(t_sim))
    V = zeros(length(t_sim))
    for i = 1:length(t_sim)
        R[i] = (sqrt(X[1, i]^2+X[2, i]^2+X[3, i]^2)-1)*Re
        V[i] = sqrt(X[8, i]^2+X[9, i]^2+X[10, i]^2)
    end
    Plots.plot(t_sim, R, label="Spacecraft Altitude", show=true, color=:red)
    Plots.plot!(twinx(), V, label="Spacecraft Velocity", show=true, color =:blue)

    xlabel!("t [s]")
    title!("Spacecraft Altitude and Velocity")
end

function plot_traj(X)
    Plots.plot(X[1, :]*Re, X[2, :]*Re, label="Spacecraft Trajectory")
    xlabel!("X [km]")
    ylabel!("Y [km]")
    title!("Spacecraft Trajectory in MCI - XY projection")
end

function plot_traj3D(X)
    using3D()
    pygui(true)
    fig = figure()
    ax = fig[:gca](projection="3d")
    Plots.plot(X[1, :]*Re, X[2, :]*Re, X[3, :]*Re, label="Spacecraft Trajectory")
    xlabel!("X [km]")
    ylabel!("Y [km]")
    #zlabel!("Z [km]")
    plot_sphere(100, Re)
    title!("Spacecraft Trajectory in MCI")
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
    Plots.plot(t_sim, R, label="Spacecraft Altitude", show=true)
    xlabel!("t [s]")
    ylabel!("altitude [km]")
    title!("Spacecraft Altitude")
end

function plot_quaternions(X, t_sim)
    Plots.plot(t_sim, X[4, :])
    Plots.plot!(t_sim, X[5, :])
    Plots.plot!(t_sim ,X[6, :])
    Plots.plot!(t_sim, X[7, :])
    xlabel!("t [s]")
    ylabel!("Quaternions")
    title!("Spacecraft Quaternions")
end

function plot_ang_vel(X, t_sim)
    Plots.plot(t_sim, X[11, :])
    Plots.plot!(t_sim, X[12, :])
    Plots.plot!(t_sim, X[13, :])
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

function plot_sphere(n, Re)
    u = range(0.0,stop = 2*pi, length = n);
    v = range(0.0,stop = pi, length= n);
    x = Re*cos.(u) .* sin.(v)'
    y = Re*sin.(u) .* sin.(v)'
    z = Re*ones(n).*cos.(v)'
    Plots.plot!(x, y, z, color = :brown, legend = false)
end

function plot_attack_angle(X, t_sim)
    Re = 3389.5
    n = length(t_sim)
    α = zeros(n)
    for i =1:1:n
        q = X[4:7, i]
        v = X[8:10, i]
        r = X[1:3, i]
        ω_mars = [0; 0; 7.095*10^(-5)]
        v_rel = (v-cross(ω_mars, r*Re)) #velocity of spacecraft wrt atm
        v_body = qrot(qconj(q), v_rel) #velocity of spacecraft wrt atm in body frame
        α[i] = acos(v_body[3]/norm(v_body)) #radians
        α[i] = α[i]*180/pi
    end
    Plots.plot(t_sim, α)
end
