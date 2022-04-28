module Geometry
export HexPoints
# plot hex
function Hex(N::Int64)
    #plot hex bx N plots in SC in 111 direction
    # todo: expand circle by circle, not 3 circles by 3 circles
    xyzs = Array{Int16}(undef,0,3)
    xyzs = ImportLatiice(xyzs,Int16(0),Int16(0),Int16(0),N)
    if size(xyzs)[1] >= N
        return xyzs
    end
    R = Int16(1)
    while true
        x = -Int16(2)*R
        y = Int16(2)*R
        z = Int16(0)
        for _ in Int16(1):R
            y = y - Int16(2)
            z = z + Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int16(1):R
            y = y - Int16(2)
            x = x + Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int16(1):R
            x = x + Int16(2)
            z = z - Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int16(1):R
            y = y + Int16(2)
            z = z - Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int16(1):R
            y = y + Int16(2)
            x = x - Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int16(1):R
            x = x - Int16(2)
            z = z + Int16(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        R += Int16(1)
    end
end


function ImportLatiice(xyzs::Matrix{Int16},x::Int16,y::Int16,z::Int16,N::Int64)
    xyzs = vcat(xyzs, [x y z])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs = vcat(xyzs, [x-Int16(1) y+Int16(1) z+Int16(1)])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs = vcat(xyzs, [x-Int16(1) y-Int16(1) z+Int16(1)])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs
end


#output xyzs in lammps dump format
function Dump(xyzs::Matrix{Int16}, filename::String, nstep::Int64)
    f = open(filename, "a")
    write(f, "ITEM: TIMESTEP\n")
    write(f, "$(nstep)\n")
    write(f, "ITEM: NUMBER OF ATOMS\n")
    write(f, "$(size(xyzs)[1]) \n")
    write(f, "ITEM: BOX BOUNDS pp pp pp\n")
    write(f, "-100 100\n")
    write(f, "-100 100\n")
    write(f, "-100 100\n")
    write(f, "ITEM: ATOMS id type x y z\n")
    for i in 1:size(xyzs)[1]
        write(f, "$(i) 1 $(xyzs[i,1]) $(xyzs[i,2]) $(xyzs[i,3]) \n")
    end
    close(f)
end


function HexPoints(N::Int64, centerCood::Vector{Int16}, direction::Vector{Int16})
    xyzs = Hex(N)
    for i in 1:3
        xyzs[:,i] *= direction[i]
        xyzs[:,i] .+= centerCood[i]
    end
    xyzs
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    filename = "/mnt/c/Users/XUKE/Desktop/hex.dump"
    f = open(filename, "w")
    close(f)
    for n = 1:100
        xyzs = Geometry.HexPoints(n, Vector{Int16}([3,3,3]), Vector{Int64}([1,1,-1]))
        Geometry.Dump(xyzs, filename, n)
    end
end



