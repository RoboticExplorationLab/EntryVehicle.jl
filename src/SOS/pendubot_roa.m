%% DYNAMICS PENDUBOT TRIAL (NEED TO BE MODIFIED, NEED REAL NON POLYNOMIAL DYNAMICS)

tspan = 0.0:0.01:2.0;
y0 = [0.5; 0.0; 0.5; 0.0];
[t_out,y_out] = ode45(@(t, y) cl_pendubot(t, y), tspan, y0);

figure()
plot(t_out, y_out(:, 1))
hold on
plot(t_out, y_out(:, 2))


%% MONTE CARLO ESTIMATION ROA FOR PENDUBOT SYSTEM

clc;clear;clear all; 

M = 100;
x1 = unifrnd(-1.7, 1.7, 1.0, M);
x3 = unifrnd(-1.7, 1.7, 1.0, M);
X = [x1;x3]
t_span = [0.0 1.0]
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-4);

Y = zeros(4, M);
for i=1:M
    y_0 = [X(1, i); 0.0; X(2, i); 0.0] %x2 and x4 set to 0 (slice of system)
    [t_out,y_out] = ode45(@(t, y) cl_pendubot(t, y), t_span, y_0)
    if abs(y_out(end, 1)) < 1e-2
        Y(:, i) = y_0;
    end
    i
end

figure()
scatter(Y(1, :), Y(3, :))
xlabel("x1")
ylabel("x3")
title("region of attraction")
xlim([-1.7, 1.7])
ylim([-1.7, 1.7])


function dy = cl_pendubot(t, y)
    % y1 and y3 angular positions of 2 links
    % actuated with torque only the first link
    dy = zeros(4, 1);
    dy(1) = y(2);
    dy(2) = 782.0*y(1)+135.0*y(2)+689.0*y(3)+90.0*y(4);
    dy(3) = y(4);
    dy(4) = 279.0*(y(1)*(y(3)^2))-1425.0*y(1)-257.0*y(2)+(273.0*(y(3)^3))-1249.0*y(3)-171.0*y(4);
end
