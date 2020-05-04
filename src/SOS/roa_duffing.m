%% ROA estimation for the DUFFING OSCILLATOR SYSTEM
% Need controller and ROA comparison

clc;clear;clear all; 
a = -1.0;
b = 1.0;
g = 0.0;
d = 10.0;
o = 1.0;
p = [a;b;g;d;o];

%% Test Plot Non linear Dynamics Behavior

tspan = 0.0:0.01:20.0;
y0 = [1.0; 1.0];
[t_out,y_out] = ode23(@(t, y) duffing(t, y, p, 0.0), tspan, y0);

figure()
plot(y_out(:, 1), y_out(:, 2))
title("phase portrait duffing")
xlabel("x1")
ylabel("x2")

%% Monte Carlo Bassin of Attraction (for the fixed found controller)

M = 1000;
x1 = unifrnd(-50.0, 50.0, 1.0, M);
x2 = unifrnd(-50.0, 50.0, 1.0, M);
X = [x1;x2];
t_span = 0.0:0.01:200.0;

Y = zeros(2, M);
for i=1:M
    y_0 = X(:, i);
    [t_out,y_out] = ode23(@(t, y) duffing(t, y, p, 0.0), t_span, y_0);
    if abs(y_out(end, 1)-1.0) < 1e-2 && abs(y_out(end, 2)) < 1e-2
        Y(:, i) = y_0;
    end
    i
end

figure()
scatter(Y(1, :), Y(2, :))
xlabel("x1")
ylabel("x2")
title("region of attraction MC estimation")

%% Duffing with control

%eq point
y1 = 5.0;
y2 = 0.0;
A = [0.0 1.0; -a-3*b*(y1)^2 -d];
B = [0.0; 1.0];

Q = eye(2);
R = 1.0;

K_lqr = lqr(A, B, Q, R);

tspan = 0.0:0.01:100.0;
y0 = [0.5; 1.0];
[t_out,y_out] = ode23(@(t, y) duffing(t, y, p, -K_lqr*y), tspan, y0);

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

%% Linearization and Initial Lyapunov Function

A = [0.0 1.0; -a-3*b -d]; %linearized at (1.0,0.0) point
Q1 = eye(2);
S = lyap(A', Q1); %solution of Lyapunov eq.

%% SOS formulation of ROA initialization

epsi = 1e-3;
d = 6;
d1 = 4;
d2 = 2;
% variables
sdpvar x1 x2 beta
f = [x2; -d*x2-a*x1-b*x1^3]; %polynomial
x = [x1;x2];
V = x'*S*x;
h = [x1 x2]*[x1; x2]; %shape factor
l1 = 1e-6*h;
l2 = 1e-6*h;
m = 0.0;

%% SOS optimization loop

for i=1:1
    u = m+20.0; %Now binary search over rho to maximize rho
    l = m-1.0; %We know that one works (we want at least this next)
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
    m = value(beta)
    V = vv{1}'*QQ{1}*vv{1}+l1;
    sdisplay(V)
end

%% Plot the different ROA estimates 

V1 = sdisplay(V);
L2=strrep(strrep(V1,'*','.*'),'^','.^');V3=cell2mat((L2));
[x1,x2]=meshgrid([-3:0.01:3],[-3:0.01:3]);
hold on 
%figure()
contour(x1,x2,eval(V3),[1 1], 'Color', 'r') %to see the level set
camlight; lighting gouraud


%% Functions

function dy = duffing(t, y, p, u)
    dy = zeros(2, 1);
    a = p(1);
    b = p(2);
    g = p(3);
    d = p(4);
    o = p(5);
    dy(1) = y(2);
    dy(2) = g*cos(o*t)-d*y(2)-a*y(1)-b*y(1)^3 + u;
end

%Additional functions used for the "shape factor" method

function [sol,v_sol,Q_sol] = step_4(beta,V,f,h,l1,l2,x,d1,d2) 
    [s1,c1] = polynomial(x,d1);
    [s2,c2] = polynomial(x,d2);
    dVdt = jacobian(V,x)*f;
    D = -(dVdt+l2)+s1*(V-1.0);
    E = (h-beta)*s2 + (1.0-V);
    F = [sos(s1),sos(s2),sos(D),sos(E)];
    ops = sdpsettings('solver','mosek','verbose',0);
    [sol,v_sol,Q_sol]=solvesos(F,[],ops,[c1;c2]);
end

function [sol,v_sol,Q_sol] = step_5(beta,s1,s2,f,h,l1,l2,x,d)
    [V,p] = polynomial(x,d);
    dVdt = jacobian(V,x)*f;
    D = -(l2+dVdt)+s1*(V-1);
    E = (h-beta)*s2+ (1-V);
    F = [sos(V-l1),sos(D),sos(E)];
    ops = sdpsettings('solver','mosek','verbose', 0.0);
    [sol,v_sol,Q_sol]=solvesos(F,-beta,ops,[beta;p]);
end