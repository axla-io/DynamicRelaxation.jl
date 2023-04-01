
function set_rotation_vels!(dr, ω, n, u_t)

     for i in 1:n
        dr_id = 4*(i-1)+1
        ω_id = 3*(i-1)+1

        dr_i = SVector{3, u_t}(@view([dr_id:dr_id+3]))
        ω_i = SVector{3, u_t}(@view([ω_id:ω_id+2]))

        # Do rotation according to Rucker, “Integrating Rotations Using Nonunit Quaternions.”
        
     end
end