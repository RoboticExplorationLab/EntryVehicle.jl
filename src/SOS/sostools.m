% This file uses SOStools for ROA problems involving SOS programming
% VAN DER POL OSCILLATOR EXAMPLE HERE
% This file enables to compare the results obtained from SOStools (here)
% to the ones obtained using YALMIP 

clc;clear;clear all; 

%% Syntax example for SOStools
syms x y;
p = 4*x^4*y^6+x^2-x*y^2+y^2;
[Q,Z]=findsos(p,"rational");

%% Test on SOStools

d = 4;
t = 2.3;
syms x1 x2;
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
S = lyap(A', eye(2));
V0 = [x1 x2]*S*x;
V = V0;
prog = sosprogram([x1;x2]);
VEC = monomials([x1; x2],0:1:d); %create vector of monomials up to degree 4
vec = monomials([x1; x2],0:1:d/2);
[prog,L] = sossosvar(prog,vec);
dVdt = [diff(V,x1), diff(V,x2)]*f;
prog = sosineq(prog,-dVdt+L*(V-t));
prog = sossolve(prog);

%% Cheng SOS paper comparison using SOStools

d = 4; %only even integers
epsi = 1e-3;
syms x1 x2;
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
Q1 = eye(2);
Q2 = [1.0 0.0; 0.0 2.0];
Q3 = [5.0 0.0; 0.0 2.0];
S = lyap(A', Q3);
V0 = [x1 x2]*S*x;
V = V0;
l1 = 1e-6*[x1 x2]*x;
l2 = l1;
m = 0.0;
for i=1:1
    u = m +10.0;
    l = m;
    while abs(u-l)>epsi
        t = (u+l)/2;
        prog = sosprogram([x1;x2]);
        vec = monomials([x1; x2],0:1:d/2);
        [prog,L] = sossosvar(prog,vec);
        dVdt = [diff(V,x1), diff(V,x2)]*f;
        prog = sosineq(prog,-(dVdt-l2)+L*(V-t)); %l2 not really necessary here
        solver_opt.solver = 'sedumi';
        solver_opt.params.tol = 1e-9;
        prog = sossolve(prog, solver_opt);
        L_sol = sosgetsol(prog,L);
        if L_sol == 0.0
            u = t;
        else 
            l = t;
        end
    end
end

%% Cheng comparison "shape factor" method using SOStools (trial 2)
%No bissection for step 2 but simply solving instead.
%WORKING

clc;clear;clear all; 

d = 6; %only even integers
d1 = 2;
d2 = 4;
epsi = 1e-3;
syms x1 x2 gam;
x = [x1;x2];
f = [-x2; x1+((x1^2)-1)*x2]; %((x1^2)-1)*x2]; %non-linear dynamics equation
A = [0.0 -1.0; 1.0 -1.0]; %linearization at point (0.0,0.0) stable equi point
Q1 = eye(2);
Q2 = [1.0 0.0; 0.0 2.0];
Q3 = [5.0 0.0; 0.0 2.0];
S = lyap(A', Q3);
%V0 = [x1 x2]*S*x;
%V = V0;
V = (0.65217*x1^2-0.43478*x1*x2+0.43478*x2^2);
l1 = 1e-6*[x1 x2]*x;
l2 = l1;
h = [x1 x2]*x;
m = 0.0;
vec1 = monomials([x1; x2],0:1:d1/2);
vec2 = monomials([x1; x2],0:1:d2/2);
vec3 = monomials([x1; x2],0:1:d/2);
for i=1:1
    u = m +10.0;
    l = m;
    while abs(u-l)>epsi
        t = (u+l)/2
        prog = sosprogram([x1;x2]);
        [prog,s1] = sossosvar(prog,vec1);
        [prog,s2] = sossosvar(prog,vec2);
        dVdt = [diff(V,x1), diff(V,x2)]*f;
        prog = sosineq(prog,-(dVdt+l2)+s2*(V-1.0)); %l2 not really necessary here
        prog = sosineq(prog, (h-t)*s1+(1.0-V));
        solver_opt.solver = 'sedumi';
        solver_opt.params.tol = 1e-12;
        prog = sossolve(prog, solver_opt);
        s1_sol = sosgetsol(prog,s1);
        s2_sol = sosgetsol(prog,s2);
        if s1_sol == 0.0 || s2_sol == 0.0
            u = t;
        else 
            l = t;
            s1_t = s1_sol;
            s2_t = s2_sol;
        end
    end
    prog1 = sosprogram([x1,x2], [gam]);
    [prog1,Vv] = sossosvar(prog1,vec3);
    dVdt = [diff(Vv,x1), diff(Vv,x2)]*f;
    K = -(dVdt+l2)+(s2_t*(Vv-1.0));
    prog1 = sosineq(prog1,((h-gam)*s1_t)+(1.0-Vv));
    prog1 = sosineq(prog1,K); 
    prog1 = sosineq(prog1, Vv-l1);
    prog1 = sossetobj(prog1, -gam);
    solver_opt.solver = 'sedumi';
    solver_opt.params.tol = 1e-9;
    prog1 = sossolve(prog1, solver_opt);
    V_sol = sosgetsol(prog1,Vv);
    V = V_sol;
    m = sosgetsol(prog1,gam)
end


%% Plot - need to be adapted to SOStools variables

V1 = sdisplay(V);
L2=strrep(strrep(V1,'*','.*'),'^','.^');V3=cell2mat((L2));
[x1,x2]=meshgrid([-3:0.01:3],[-3:0.01:3]);
hold on 
%figure() %remove figure() to plot on top of the other methods regions
%estimations
contour(x1,x2,eval(V3),[1 1], 'Color', 'r') %to see the level set
camlight; lighting gouraud
