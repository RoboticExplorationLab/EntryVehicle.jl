# Entry Vehicle Geometry

#TO DO: add more missions geometry

abstract type Geometry{T} end

"""
Geometry
Holds all information required to uniquely describe an entry vehicle geometry.
So far using only Spherecone geometry
"""

struct SphereConeGeometry{T} <: Geometry{T}
    δ::T # semi-apex angle cone part in radians (rad)
    r_cone::T # radius base cone part in meters (M=m)
    r_min::T # radius sphere part in meters (M=m)
end

function get_r_sphere(geometry::SphereConeGeometry{T}) where {T}
    return geometry.r_min/(cos(geometry.δ))
end

function get_A_ref(geometry::SphereConeGeometry{T}) where {T}
    # Return Reference Area of Sphere Cone
    return π*geometry.r_cone^2
end

function get_L_ref(geometry::SphereConeGeometry{T}) where {T}
    # Return Reference Length of Sphere Cone
    return geometry.r_cone
end

function PhoenixSphereConeGeometry()
    #Phoenix geometry from
    #Aerodynamics for the Mars Phoenix Entry Capsule
    #by Karl T. Edquist†, Prasun N. Desai*, and Mark Schoenenberger‡ - page 3
    SphereConeGeometry{Float64}(70.0*pi/180, 1.325, 0.2)
end

function MSLSphereConeGeometry()
    #Find Numbers on that one
    SphereConeGeometry{Float64}()
end


# Space Entry Vehicle Sturcture

abstract type EntryVehicle{T} end

#Might make distinctions between Capsule and Shuttle in the future

"""
Vehicle
Holds all information required to uniquely describe a Vehicle. A Vehicle is
described by its geometry (type Geometry{T}) and by its mass properties
"""

struct Vehicle{T} <: EntryVehicle{T}
    geometry::Geometry{T}     # Geometry of the Vehicle
    m::T                      # total mass of the vehicle - kg
    J::Matrix{T}              # inertia matrix of the vehicle
    Jinv::Matrix{T}           # inverse inertia matrix of the vehicle
    r_com::Vector{T}          # position COM wrt base cone in meters
end

function PhoenixVehicle()
    PhoenixGeometry = PhoenixSphereConeGeometry()
    m = 600.0  # Mass kg
    J = [184.180 0.0 0.0;
         0.0 184.180 0.0;
          0.0 0.0 230.0]
    J_inv = [0.00542947 0.0 0.0;
             0.0 0.00542947 0.0;
             0.0 0.0 0.00434783]
    r_com = [0.001,0.,-0.189]  # COM position in meters
    return Vehicle{Float64}(PhoenixGeometry, m, J, J_inv, r_com)
end
