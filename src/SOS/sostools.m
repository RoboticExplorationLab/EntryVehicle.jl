% This file uses SOStools for ROA problems involving SOS programming

clc;clear;clear all; 

%% Syntax example
syms x y;
p = 4*x^4*y^6+x^2-x*y^2+y^2;
[Q,Z]=findsos(p,"rational");

%% Test

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
VEC = monomials([x1; x2],1:1:d); %create vector of monomials up to degree 4
[prog,L] = sossosvar(prog,VEC);
dVdt = [diff(V,x1), diff(V,x2)]*f;
prog = sosineq(prog,-dVdt+L*(V-t));
prog = sossolve(prog);

