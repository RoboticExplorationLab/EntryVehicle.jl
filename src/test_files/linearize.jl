using ForwardDiff
using LinearAlgebra
using StaticArrays

include("quaternions.jl")

a=1

function f(u)
    I1 = 100.0 #p[1] #assume principal axes
    I2 = 110.0 #p[2]
    I3 = 300.0 #p[3]
    #τ = F(t) #see how to deal with that
    τ = [0.0;0.0;0.0]
    #@show(τ)
    I = Diagonal([I1; I2; I3])
    I_inv = Diagonal([1/I1; 1/I2; 1/I3])
    #controller
    q_ref = [0.0; 1.0; 0.0; 0.0]
    kd = 30.0 #30.0
    kp = -20.0 #20.0
    q_err = qmult(qconj(u[1:4]), q_ref) #perfect measurements
    #@show(qconj(u[1:4]))
    τ_c = -kd*u[5:7]-kp*q_err[2:4]
    #τ_c = [0.0;0.0;0.0]
    #update
    u[5:7] = I_inv*(τ-cross(u[5:7], I*u[5:7])+τ_c)
    u[1:4] = 0.5*qmult(u[1:4], [0.0; u[5:7]])
    return u
end

f(u)

ForwardDiff.jacobian(u -> f(u[1], u[2], u[3], u[4], u[5], u[6], u[7]), u); # g = ∇f

u = [1.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.5]
g(u)


f(x,y)=[x^2+y^3-1,x^4 - y^4 + x*y]
ForwardDiff.jacobian(x ->f(x[1],x[2]), [1.0, 1.0])

f(x::Vector) = sum(sin, x) + prod(tan, x) * sum(sqrt, x);
g = x -> ForwardDiff.gradient(f, x); # g = ∇f
x = rand(5)
g(x)


I1 = 100.0
I2 = 110.0
I3 = 300.0
I = Matrix(Diagonal([I1; I2; I3]))
I_inv = Matrix(Diagonal([1/I1; 1/I2; 1/I3]))
u = [1.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.5]
τ = [1.0; 0.0; 0.0]
f(u) = [0.5*qmult(u[1:4], [0.0; u[5:7]]); I_inv*(τ-cross(u[5:7], I*u[5:7]))]

a = [1.0; 0.0; 0.0; 0.0; 0.1; 0.1; 0.5]
A = ForwardDiff.jacobian(u ->f(u), a)
