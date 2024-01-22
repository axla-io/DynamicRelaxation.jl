struct CoordinateSystem
    x::SVector{3,Float64}
    y::SVector{3,Float64}
    z::SVector{3,Float64}
end

function CoordinateSystem()
    CoordinateSystem(SVector{3,Float64}(1, 0, 0), SVector{3,Float64}(0, 1, 0), SVector{3,Float64}(0, 0, 1))
end

# Struct for storing element properties
abstract type AbstractElement end
struct Beam<:AbstractElement
    E::Float64    # [Pa]
    A::Float64    # [m^2]
    Iy::Float64   # [m^4]
    Iz::Float64   # [m^4]
    G::Float64    # [Pa]
    It::Float64   # [m^4]
    l_init::Float64   # [m]
    cs::CoordinateSystem
end

struct Bar<:AbstractElement
    E::Float64    # [Pa]
    A::Float64    # [m^2]
    l_init::Float64   # [m]
    cable::Bool 
end

function ElementProperties(E, A, Iy, Iz, G, It, l_init)
    Beam(E, A, Iy, Iz, G, It, l_init, CoordinateSystem())
end

function ElementProperties(E, A, l_init, cable)
    Bar(E, A, l_init, cable)
end