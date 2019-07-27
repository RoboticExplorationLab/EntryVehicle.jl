#tests plots in Python
using PyPlot
using3D()
pygui(true)

fig = figure()
ax = fig[:gca](projection="3d")

#N = 10
#x,y,z,u,v,w = [randn(N) for _ in 1:6]
#ax[:quiver](x,y,z, u,v,w)
#xlim([minimum(x), maximum(u+x)])
#ylim([minimum(y), maximum(v+y)])
#zlim([minimum(z), maximum(w+z)])

n = 100
u = range(0.0,2*π,length=n);
v = range(0.0,π,length=n);

Re = 6378.1
x = Re*cos.(u) * sin.(v)';
y = Re*sin.(u) * sin.(v)';
z = Re*ones(n) * cos.(v)';

surf(x,y,z)
