
function set_rotation_vels!(dr, ω, n, u_t)

     for i in 1:n
        dr_id = 4*(i-1)+1
        ω_id = 3*(i-1)+1

        dr_i = SVector{3, u_t}(@view([dr_id:dr_id+3]))
        ω_i = SVector{3, u_t}(@view([ω_id:ω_id+2]))

        # Do rotation according to Rucker, “Integrating Rotations Using Nonunit Quaternions.”
        
     end
end

function scross(a::SVector{3, T}, b::SVector{3, T}) where T<:Real
    a1, a2, a3 = a
    b1, b2, b3 = b

    SVector{3, T}(a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1)
end