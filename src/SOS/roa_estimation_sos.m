% Region of Attraction Inverted Pendulum

%%

clc;clear;clear all; 


%% Test Plot Non linear Dynamics Behavior

tspan = 0.0:0.01:20.0;
y0 = [3*pi/4; 0.0];
[t_out,y_out] = ode45(@(t, y) inv_pend(t, y, 0.0), tspan, y0);

figure()
plot(t_out, y_out(:, 1))

%% Linearized Dynamics about theta = 0.0 (upright pendulum)

g=10.0;
l=1.0;
m=2.0;
b = 2.0;
A = [0.0 1.0; g/l -b/(m*l^2)];
B = [0.0; 1/(m*l^2)];

Q = 0.1*[1.0 0.0; 0.0 0.1];
R = 100.0;

K_lqr = lqr(A, B, Q, R);


tspan = 0.0:0.001:10.0;
y0 = [0.88; 0.0];
[t_out,y_out] = ode45(@(t, y) inv_pend(t, y, -K_lqr*y), tspan, y0);

u = zeros(size(y_out(:,1)));
for i=1:length(y_out(:, 1))
    u(i) = -K_lqr*y_out(i, :)';
end

figure()
plot(t_out, y_out(:, 1))

figure()
plot(t_out, u);

%% Monte Carlo Bassin of Attraction (for the fixed found controller)

M = 10000;
x1 = unifrnd(0.0, 2*pi, 1.0, M);
x2 = unifrnd(-5.0, 5.0, 1.0, M);
X = [x1;x2];
t_span = 0.0:0.01:5.0;

Y = zeros(2, M);
for i=1:M
    y_0 = X(:, i);
    [t_out,y_out] = ode45(@(t, y) inv_pend(t, y, -K_lqr*y), tspan, y_0);
    if y_out(end, 1) < 1e-3 
        Y(:, i) = y_0;
    end
end

figure()
scatter(Y(1, :), Y(2, :))
xlabel("x1")
ylabel("x2")
title("region of attraction")

%% ROA estimation using SOS 

% variables
sdpvar x1 x2 rho
% nonlinear system dx(t)/dt=f(x)
%f=[-x1+(1+x1)*x2;-(1+x1)*x1];
%f=inv_pend(0.0, [x1;x2], -K_lqr*[x1;x2]); %time value does not matter : autonomous sys
g = 10.0;
l = 1.0;
b = 2.0;
m = 2.0;
f = [x2; (g/l)*sin(x1)-K_lqr*[x1;x2]/(m*l^2)-(b/(m*l^2))*x2];
x = [x1;x2];
% Polynomial Lyapunov Function of order 4 
[V,c] = polynomial(x,4);
[L,p] = polynomial(x,4); 
% dV(x)/dt
dVdt = jacobian(V,x)*f;
D = -dVdt+L*(V-rho);
% SOS COnditions
F = [sos(V),sos(D),sos(L),sum(c)==1.0,c(1)==0.0];

%% 

%SDP solver
ops = sdpsettings('solver','mosek');
% SOlve SOS program
[sol,v,Q]=solvesos(F,-rho,ops,[c;p]);

%% Functions 

function dy = inv_pend(t, y, u)
    dy = [0.0; 0.0];
    %if u >= 5.0 || u <= -5.0
    %    dy(1) = 0.0;
    %    dy(2) = 0.0;
    %else
    m = 2.0;
    g = 10.0;
    l = 1.0;
    b = 2.0;
    dy(1) = y(2);
    dy(2) = (g/l)*sin(y(1))+u/(m*l^2)-b/(m*l^2)*y(2);
    %end
end