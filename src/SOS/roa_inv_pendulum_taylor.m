% Region of Attraction Inverted Pendulum Example with Taylor expansion

%%

clc;clear;clear all; 

%% Test Plot Non linear Dynamics Behavior

tspan = 0.0:0.01:50.0;
y0 = [3*pi/4; 0.0];
[t_out,y_out] = ode45(@(t, y) inv_pend(t, y, 0.0), tspan, y0);

figure()
plot(t_out, y_out(:, 1))

%% Linearized Dynamics about theta = 0.0 (upright pendulum) unstable 
% point.

g=1.0;
l=1.0;
m=1.0;
b = 0.5;
A = [0.0 1.0; g/l -b/(m*l^2)];
B = [0.0; 1/(m*l^2)];

Q = 1e-6*eye(2);
R = 1e8;

K_lqr = lqr(A, B, Q, R);
%K_lqr = [100.0 67.0];

tspan = 0.0:0.01:100.0;
y0 = [pi/2; -500.0];
[t_out,y_out] = ode23(@(t, y) inv_pend(t, y, -K_lqr*y), tspan, y0);

u = zeros(size(y_out(:,1)));
for i=1:length(y_out(:, 1))
    u(i) = -K_lqr*y_out(i, :)';
end

figure()
plot(t_out, y_out(:, 1))
hold on
plot(t_out, y_out(:, 2))

figure()
plot(t_out, u);

%% Monte Carlo Bassin of Attraction (for the fixed found controller)

M = 1000;
x1 = unifrnd(0.0, 2*pi, 1.0, M);
x2 = unifrnd(-50.0, 50.0, 1.0, M);
X = [x1;x2];
t_span = 0.0:0.01:200.0;

Y = zeros(2, M);
for i=1:M
    y_0 = X(:, i);
    [t_out,y_out] = ode23(@(t, y) inv_pend(t, y, -K_lqr*y), tspan, y_0);
    if abs(y_out(end, 1)) < 1e-3 && abs(y_out(end, 2)) < 1e-3
        Y(:, i) = y_0;
    end
    i
end

figure()
scatter(Y(1, :), Y(2, :))
xlabel("x1")
ylabel("x2")
title("region of attraction")

%% ROA estimation using SOS with Taylor expansion of pendulum

clc;clear;clear all; 

epsi = 1e-3;
d = 6;
d1 = 4;
d2 = 2;
% variables
sdpvar x1 x2 beta
s = x1-((x1^3)/6)+((x1^5)/5);
%f = [x2; (g/l)*s-K_lqr*[x1;x2]/(m*l^2)-(b/(m*l^2))*x2]; %polynomial
f = [x2; (g/l)*s-(b/(m*l^2))*x2]; %polynomial with no control
x = [x1;x2];
Q1 = eye(2);
Q2 = [1.0 0.0; 0.0 2.0];
Q3 = [5.0 0.0; 0.0 2.0];
S = lyap(A', Q1);
V = x'*S*x;
h = [x1 x2]*[x1; x2]; %shape factor
l1 = 1e-6*h;
l2 = 1e-6*h;
va = 0.0;

%% ROA of LQR controller using SOS


for i=1:1
    u = va+20.0; %Now binary search over rho to maximize rho
    l = va-1.0; %We know that one works (we want at least this next)
    while abs(u-l)>1e-3
        t = (u+l)/2
        [s,v,Q] = step_4(t,V,f,h,l1,l2,x,d1,d2);
        result = s.problem
        if result == 0
            l = t;
            s1 = v{1}'*Q{1}*v{1};
            s2 = v{2}'*Q{2}*v{2};
        else
            u = t;
        end
    end
    %s1 = clean(s1/(sum(coefficients(s1,x))), 1e-9);
    %s2 = clean(s2/(sum(coefficients(s2,x))), 1e-9);
    [ss,vv,QQ] = step_5(beta,s1,s2,f,h,l1,l2,x,d);
    ss.problem
    va = value(beta)
    V = vv{1}'*QQ{1}*vv{1}+l1;
    sdisplay(V)
end



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
    m = 1.0;
    g = 1.0;
    l = 1.0;
    b = 0.5;
    dy(1) = y(2);
    dy(2) = (g/l)*sin(y(1))+u/(m*l^2)-b/(m*l^2)*y(2);
    %end
end