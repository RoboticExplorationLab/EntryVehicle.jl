import Base: copy, reset

using StaticArrays

abstract type Geometry end

"""
Geometry
Holds all information required to uniquely describe an entry vehicle geometry.
So far using only Spherecone geometry
"""

struct SphereconeGeometry <: Geometry
    δ::Float64 # semi-apex angle cone part in radians
    r_cone::Float64 # radius base cone part in meters
    r_min::Float64 # radius sphere part in meters
end

function r_sphere(SphereconeGeometry)
    return SphereconeGeometry.r_min/(cos(SphereconeGeometry.δ))
end

function A_ref(SphereconeGeometry)
    return pi*SphereconeGeometry.r_cone^2
end

function L_ref(SphereconeGeometry)
    return SphereconeGeometry.r_cone
end

################################################################################
################################################################################
"""
Vehicle
Holds all information required to uniquely describe a Vehicle. A Vehicle is
described by its geometry (type Geometry) and by its mass properties
"""

struct Vehicle
    name_vehicle::String #name vehicle
    geometry::Geometry
    m::Float64 # mass of the vehicle - kg
    J::SMatrix{3,3} # inertia matrix of the vehicle
    Jinv::SMatrix{3,3} # inertia matrix inverse of the vehicle
    r_com::SVector{3} # position COM wrt base cone in meters
end


geo_phoenix=SphereconeGeometry(70.0*pi/180, 1.325, 0.2)
phoenix=Vehicle("phoenix",geo_phoenix,600.0,SMatrix{3,3}([184.180 0.0 0.0;
                                        0.0 184.180 0.0;
                                        0.0 0.0 230.0]),
                    SMatrix{3,3}([0.00542947 0.0 0.0;
                                        0.0 0.00542947 0.0;
                                        0.0 0.0 0.00434783]),
                            SVector{3}([0.001,0.,-0.189]))

################################################################################
################################################################################

"""
Environment
Holds all information required to uniquely describe a planetary environment.
Earth and Mars will be primarily used
"""


struct Environment
    name_env::String # name of the environment
    R_p::Float64 # Radius of the central planet - m
    μ_p::Float64 # Gravitational Parameter for central planet - m3.s-2
    M_p::Float64 # Mass central planet - kg
    J2::Float64 # J2 parameter value
    ρ0::Float64 # Atmospheric Density parameter (exp model) - kg.m-3
    H::Float64 # Characteristic height in exp model - m
end

Earth=Environment("earth",6378.135*1e3,3.986004418*1e14,5.972e24,1082.48e-6,
                    1.22,8.0e3);

Mars=Environment("mars",3389.5*1e3,4.282837e13,6.39*1e23,1.96e-3,0.0158,
                9.35458*1e3);
