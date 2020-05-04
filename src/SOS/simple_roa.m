% Test on simple system. ROA known

clc;clear;clear all; 

epsi = 1e-3;
d = 6;
sdpvar x1 x2 rho
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
S = lyap(A', eye(2));
V0 = x'*S*x;
V = V0;
imp = 10.0; %improvement from one outer iteration to the other
m = 0.0;

%% %while abs(imp) > 1e-1 %no need for abs because guaranteed to improve at each step rho
for i=1:4
    u = m+20.0; %Now binary search over rho to maximize rho
    l = m-20.0; %We know that one works (we want at least this next)
    while abs(u-l)>epsi
        t = (u+l)/2
        [s3,v3,Q3] = step_2(t,V,f,x,d);
        result = s3.problem
        if result == 0
            l = t;
        else
            u = t;
        end
    end
    L = v3{1}'*Q3{1}*v3{1};
    L = clean(L/(sum(coefficients(L,x))), 1e-5);
    sdisplay(L)
    [s1,v1,Q1] = step_1(rho,L,f,x,d);
    res = s1.problem
    if res == 1 || res == 4
        temp = t;
        break
    end
    V = clean(v1{1}'*Q1{1}*v1{1}, 1e-5);
    temp = value(rho)
    imp = temp - m
    m = temp;
    %l
    %u
end

%% Plot result simple system (3)

V1 = sdisplay(V);

L2=strrep(strrep(V1,'*','.*'),'^','.^');V3=cell2mat((L2));

[x1,x2]=meshgrid([-3:0.01:3],[-3:0.01:3]);
%surf(x1,x2,eval(V3),'FaceColor','red','FaceAlpha',0.2,'EdgeColor','none','FaceLighting','phong');hold on;grid on;
hold on 
%figure()
contour(x1,x2,eval(V3),[temp temp]) %to see the level set
%xlim([-3.0 3.0]);
%ylim([-3.0, 3.0]);
camlight; lighting gouraud

%% Step by step solving
clc;clear;clear all; 

epsi = 1e-2;
d = 4;
sdpvar x1 x2 rho
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
S = lyap(A', eye(2));
V0 = x'*S*x;
V = V0;
imp = 10.0; %improvement from one outer iteration to the other
m = 0.0;
g = [x1+1;1-x1;x2+1;1-x2];

%%
t = 2.34;
[b,c] = polynomial(x,d);
[s1,c1] = polynomial(x,d); [s2,c2] = polynomial(x,d); [s3,c3] = polynomial(x,d); [s4,c4] = polynomial(x,d);
[r1,d1] = polynomial(x,d); [r2,d2] = polynomial(x,d); [r3,d3] = polynomial(x,d); [r4,d4] = polynomial(x,d);
dVdt = jacobian(V,x)*f;
D = -dVdt+b*(V-t);
sdisplay(dVdt);
sdisplay(D);
F = [sos(b-[s1 s2 s3 s4]*g),sos(D-[r1 r2 r3 r4]*g),sos(s1),sos(s2),sos(s3),sos(s4),sos(r1),sos(r2),sos(r3),sos(r4)];
ops = sdpsettings('solver','mosek','verbose',1);
[sol,v_sol,Q_sol]=solvesos(F,[],ops,[c,c1,c2,c3,c4,d1,d2,d3,d4]);
sol



%[s3,v3,Q3] = step_3(t,V,g,f,x,d);
%[s2,v2,Q2] = step_2(t,V,f,x,d);
%s3
%s2
%L = v3{1}'*Q3{1}*v3{1} + [v3{3}'*Q3{3}*v3{3} v3{4}'*Q3{4}*v3{4} v3{5}'*Q3{5}*v3{5} v3{6}'*Q3{6}*v3{6}]*g;
%L2 = v2{1}'*Q2{1}*v2{1};
%sdisplay(L)
%sdisplay(L2)

%% 

[s1,v1,Q1] = step_1(rho,L,f,x,d,t);
s1
V = v1{1}'*Q1{1}*v1{1};
sum(coefficients(V, x))
value(rho)

%% MONTE CARLO of that

M = 10000;
x1 = unifrnd(-3.0, 3.0, 1.0, M);
x2 = unifrnd(-3.0, 3.0, 1.0, M);
X = [x1;x2];
t_span = 0.0:0.01:20.0;

Y = zeros(2, M);
for i=1:M
    y_0 = X(:, i);
    [t_out,y_out] = ode45(@(t, y) vdp(t, y), t_span, y_0);
    if abs(y_out(end, 1)) < 1e-2
        Y(:, i) = y_0;
    end
end

figure()
scatter(Y(1, :), Y(2, :))
xlabel("x1")
ylabel("x2")
title("region of attraction")
xlim([-3.0, 3.0])
ylim([-3.0, 3.0])

%% Trying with L first (in this case, it seems like we are guaranteed to have rho bigger each time

clc;clear;clear all; 

epsi = 1e-3;
d = 2;
sdpvar x1 x2 rho
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
%S = lyap(A', eye(2));
L0 = 2210.17476413-66.7035713986*x1-13.7245988961*x2+653.524082004*x1^2-36.0637468352*x1^2*x2+817.262138418*x2^2-173.027905548*x1*x2-32.8920633491*x1^3-37.7990459176*x1*x2^2-48.6844111857*x2^3+1104.32552518*x1^4-68.4605236493*x1^3*x2+821.941204664*x1^2*x2^2-210.236919688*x1*x2^3+1860.19422608*x2^4;
L = 1.27432484961e-05*x2+0.230659249609*x1^2+0.0613676949508*x1^2*x2+0.04930977377*x2^2-0.140240241283*x1*x2-0.0196010454086*x1^3+0.137911053615*x1*x2^2-0.162343861625*x2^3+1.41543136901*x1^4-0.979901429523*x1^3*x2+0.701308840379*x1^2*x2^2-0.731229046559*x1*x2^3+0.437309363289*x2^4;
%L = 1.27432484961e-05*x2+0.230659249609*x1^2+0.0613676949508*x1^2*x2+0.04930977377*x2^2-0.140240241283*x1*x2;
P = rand(3,3);
%L = [1; x1; x2]'*P'*P*[1; x1; x2];
%L = clean(L/sum(coefficients(L, x)), 1e-5);
imp = 10.0; %improvement from one outer iteration to the other
m = 0.0;
level = [];
V_list = [];

%% %while abs(imp) > 1e-1 %no need for abs because guaranteed to improve at each step rho
for i=1:5
    sdisplay(L);
    [s1,v1,Q1] = step_1(rho,L,f,x,d);
    res = s1.problem
    if res == 1 || res == 4
        break
    end
    V = clean(v1{1}'*Q1{1}*v1{1}, 1e-5);
    V_list = [V_list V];
    m = value(rho)
    level = [level m];
    u = m+200.0; %Now binary search over rho to maximize rho
    l = m; %We know that one works (we want at least this next)
    while abs(u-l)>epsi
        t = (u+l)/2
        %level = [level t];
        [s3,v3,Q3] = step_2(t,V,f,x,d);
        result = s3.problem
        %if result == 4 
        %   s3
        %    sdisplay(L)
        %    sdisplay(V)
        %    m
        %    t
        %    break
        %end
        if result == 0
            l = t;
        else
            u = t;
        end
        end
    if result == 4 
        break
    end
    temp = t;
    imp = temp-m
    L = v3{1}'*Q3{1}*v3{1};
    L = clean(L/(sum(coefficients(L,x))), 1e-5);
    
end


%% functions

function dy = vdp(t, y)
    dy = zeros(2, 1);
    dy(1) = -y(2);
    dy(2) = y(1)+((y(1)^2)-1)*y(2);
end


function [sol,v_sol,Q_sol] = step_1(rho,L,f,x,d)
    [a,p] = polynomial(x,d);
    dVdt = jacobian(a,x)*f;
    D = -dVdt+L*(a-rho);
    F = [sos(a),sos(D),sum(p)==1.0, p(1)==0.0];
    ops = sdpsettings('solver','mosek', 'mosek.MSK_DPAR_ANA_SOL_INFEAS_TOL', 1.0e-10,... 
    'verbose', 0.0);
    [sol,v_sol,Q_sol]=solvesos(F,-rho,ops,[p;rho]);
end

function [sol,v_sol,Q_sol] = step_2(rho,V,f,x,d) 
    [L,c] = polynomial(x,d);
    dVdt = jacobian(V,x)*f;
    D = -dVdt+L*(V-rho);
    F = [sos(L),sos(D)];
    ops = sdpsettings('solver','mosek','verbose',0);
    [sol,v_sol,Q_sol]=solvesos(F,[],ops,c);
end

function [sol,v_sol,Q_sol] = step_3(t,V,g,f,x,d)
    [b,c] = polynomial(x,d);
    [s1,c1] = polynomial(x,d); [s2,c2] = polynomial(x,d); [s3,c3] = polynomial(x,d); [s4,c4] = polynomial(x,d);
    [r1,d1] = polynomial(x,d); [r2,d2] = polynomial(x,d); [r3,d3] = polynomial(x,d); [r4,d4] = polynomial(x,d);
    dVdt = jacobian(V,x)*f;
    D = -dVdt+b*(V-t);
    F = [sos(b-[s1 s2 s3 s4]*g),sos(D-[r1 r2 r3 r4]*g),sos(s1),sos(s2),sos(s3),sos(s4),sos(r1),sos(r2),sos(r3),sos(r4)];
    ops = sdpsettings('solver','mosek','verbose',1);
    [sol,v_sol,Q_sol]=solvesos(F,[],ops,[c,c1,c2,c3,c4,d1,d2,d3,d4]);
end
