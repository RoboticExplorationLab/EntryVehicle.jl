#Units are L=km M=kg T=hours

# TO DO : Fix discrepancies in units for different structures

"""
GravityFields
Holds all information required to uniquely describe a field of gravity
"""

# Taken From EntryGuidance.jl


abstract type AbstractGravityField{T} end

struct SphericalGravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter
end

struct J2GravityField{T} <: AbstractGravityField{T}
    μ::T #standard gravitational parameter (L^3/T^2)
    R::T #planet radius (L)
    J2::T #J2 coefficient (dimensionless)
end

function EarthGravity()
    SphericalGravityField{Float64}(398600.4418*(3600.0^2))
end

function MarsGravity()
    SphericalGravityField{Float64}(42828.37*(3600.0^2))
end


function EarthGravityJ2()
    J2GravityField{FLoat64}(398600.4415*(3600.0^2), 6378.1363, 0.1082635854e-2)
end

function MarsGravityJ2()
    J2GravityField{Float64}(42828.37*(3600.0^2), 3396.203986, 0.196045e-2)
end

function gravitational_acceleration(r::AbstractVector{T}, g::SphericalGravityField{T}) where {T}
    R = norm(r)
    a = (-g.μ/(R*R*R)).*r
end

function gravitational_acceleration(r::AbstractVector{T}, g::J2GravityField{T}) where {T}
    R = norm(r)
    R3 = R^3
    R7 = R^7
    x = r[1]
    y = r[2]
    z = r[3]
    x2y2 = x*x + y*y
    z2 = z*z
    J2 = g.J2*(g.μ*g.R*g.R) #convert from dimensionless J2 to L^5/T^2 units
    a = Diagonal([-g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(6.0*z2 - 1.5*x2y2)/R7, -g.μ/R3 + J2*(3.0*z2 - 4.5*x2y2)/R7])*r
end

"""
Environment
Holds all information required to uniquely describe a planetary environment.
Earth and Mars will be primarily used.
Environment uses GravityFields and Atmospheres and add planetary information
"""


struct Environment{T}
    R::T # Radius of the central planet - m
    μ::T # Gravitational Parameter for central planet - m3.s-2
    M::T # Mass central planet - kg
    Ω::T # Angular Velocity planet - rad.s-1
    gravity::AbstractGravityField{T} # Gravity model
    atmosphere::AbstractAtmosphere{T} # Atmospheric Model
end


function EarthEnvironment()
    Earth_Gravity = EarthGravity()
    EarthAtmosphere = EarthExponentialAtmosphere()
    env = Environment{Float64}(6378.135*1e3,3.986004418*1e14,5.972e24,7.292e-5,Earth_Gravity, EarthAtmosphere)
    return env
end

function MarsEnvironment()
    Mars_Gravity = MarsGravity()
    MarsAtmosphere = MarsExponentialAtmosphere()
    env = Environment{Float64}(3389.5*1e3,4.282837e13,6.39*1e23,7.095e-5,Mars_Gravity, MarsAtmosphere)
    return env
end


function EarthEnvironmentJ2()
    env = Environment{Float64}(6378.135*1e3,3.986004418*1e14,5.972e24,7.292e-5,EarthGravityJ2(), EarthExponentialAtmosphere())
end

function MarsEnvironmentJ2()
    env = Environment{Float64}(3389.5*1e3,4.282837e13,6.39*1e23,7.095e-5,MarsGravityJ2(), MarsExponentialAtmosphere())
end

function atmospheric_density(r::AbstractVector{T}, env::Environment{T}) where {T}
    atmospheric_density(r, env.atmosphere)
end

function atmospheric_density(r::T, env::Environment{T}) where {T}
    atmospheric_density(r, env.atmosphere)
end

function gravitational_acceleration(r::AbstractVector{T}, env::Environment{T}) where {T}
    gravitational_acceleration(r, env.gravity)
end
