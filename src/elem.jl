struct CoordinateSystem
    x::SVector{3,Float64}
    y::SVector{3,Float64}
    z::SVector{3,Float64}
end

function CoordinateSystem()
    CoordinateSystem(SVector{3,Float64}(1, 0, 0), SVector{3,Float64}(0, 1, 0), SVector{3,Float64}(0, 0, 1))
end

# Struct for storing element properties
struct ElementProperties
    E::Float64    # [Pa]
    A::Float64    # [m^2]
    Iy::Float64   # [m^4]
    Iz::Float64   # [m^4]
    G::Float64    # [Pa]
    It::Float64   # [m^4]
    l_init::Float64   # [m]
    cs::CoordinateSystem
end

function ElementProperties(E, A, Iy, Iz, G, It, l_init)
    ElementProperties(E, A, Iy, Iz, G, It, l_init, CoordinateSystem())
end
