% Code for simple inverted pendulum in polynomial state form
% IN PROGRESS

clc;clear;clear all; 
p = [2.0; 10.0; 1.0; 2.0]; %m,g,l,b

%% Test Dynamics 

tspan = 0.0:0.01:20.0;
y0 = [0.0; 1.0; 0.0; 0.0];
[t_out,y_out] = ode45(@(t, y) inv_pend_poly2(t, y, [0.0;0.0],p), tspan, y0);

figure()
plot(t_out, y_out(:, 4))

%% Linearizing 
m = p(1);
g = p(2);
l = p(3);
b = p(4);
A = [0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 (-b/m) 0.0; 0.0 0.0 0.0 (b/m)];
J = [0.0 0.0 1.0 0.0];


%% Functions Dynamics

function dy = inv_pend(t, y, u)
    dy = [0.0; 0.0];
    %if u >= 5.0 || u <= -5.0
    %    dy(1) = 0.0;
    %    dy(2) = 0.0;
    %else
    m = 10.0;
    g = 10.0;
    l = 1.0;
    b = 1.0;
    dy(1) = y(2);
    dy(2) = (g/l)*sin(y(1))+u/(m*l^2)-b/(m*l^2)*y(2);
    %end
end

function dy = inv_pend_poly(t, y, u, p)
    dy = [0.0; 0.0];
    m = p(1);
    g = p(2);
    l = p(3);
    b = p(4);
    dy(1) = (-b/m)*y(1)+u(1)/m;
    dy(2) = (-b/m)*y(2)+g+u(2)/m;
end

function dy = inv_pend_poly2(t, y, u, p)
    dy = [0.0; 0.0; 0.0; 0.0];
    m = p(1);
    g = p(2);
    l = p(3);
    b = p(4);
    dy(3) = (-b/m)*y(3)+u(1)/m;
    dy(4) = (-b/m)*y(4)+g+u(2)/m;
    dy(1) = y(3);
    dy(2) = y(4);
end