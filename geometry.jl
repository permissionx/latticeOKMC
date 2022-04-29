module Geometry
using StatsBase
export HexPoints
# plot hex

function CoordInBCC(coord::Vector{Float64})
    # To do: should shift to the nearest lattice point, but does not matter seriously
    newCoord = Vector{Int32}(undef, 3)
    newCoord[3] = Int32(round(coord[3]))
    if iseven(newCoord[3])
        for i in 1:2
            x = floor(coord[i])
            if iseven(x)
                newCoord[i] = Int32(x)
            else
                newCoord[i] = Int32(x+1)
            end
        end
    else
        for i in 1:2
            x = floor(coord[i])
            if isodd(x)
                newCoord[i] = Int32(x)
            else
                newCoord[i] = Int32(x+1)
            end
        end
    end
    newCoord
end

function Hex(N::Int64)
    #plot hex bx N plots in SC in 111 direction
    # todo: expand circle by circle, not 3 circles by 3 circles
    xyzs = Array{Int32}(undef,0,3)
    xyzs = ImportLatiice(xyzs,Int32(0),Int32(0),Int32(0),N)
    if size(xyzs)[1] >= N
        return xyzs
    end
    R = Int32(1)
    while true
        x = -Int32(2)*R
        y = Int32(2)*R
        z = Int32(0)
        for _ in Int32(1):R
            y = y - Int32(2)
            z = z + Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int32(1):R
            y = y - Int32(2)
            x = x + Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int32(1):R
            x = x + Int32(2)
            z = z - Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int32(1):R
            y = y + Int32(2)
            z = z - Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int32(1):R
            y = y + Int32(2)
            x = x - Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        for _ in Int32(1):R
            x = x - Int32(2)
            z = z + Int32(2)
            xyzs = ImportLatiice(xyzs,x,y,z,N)
            if size(xyzs)[1] >= N
                return xyzs
            end
        end
        R += Int32(1)
    end
end


function ImportLatiice(xyzs::Matrix{Int32},x::Int32,y::Int32,z::Int32,N::Int64)
    xyzs = vcat(xyzs, [x y z])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs = vcat(xyzs, [x-Int32(1) y+Int32(1) z+Int32(1)])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs = vcat(xyzs, [x-Int32(1) y-Int32(1) z+Int32(1)])
    if size(xyzs)[1] >= N
        return xyzs
    end
    xyzs
end


#output xyzs in lammps dump format
function Dump(xyzs::Matrix{Int32}, filename::String, nstep::Int64)
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


function HexPoints(N::Int64, centerCood::Vector{Float64}, direction::Vector{Int32})
    xyzs = Hex(N)
    for i in 1:3
        xyzs[:,i] *= direction[i]
    end
    nowCenter = [mean(xyzs[:,i]) for i in 1:3]
    deltaDisplace = CoordInBCC(centerCood - nowCenter)
    for i in 1:3
        xyzs[:,i] .+= deltaDisplace[i]
    end
    xyzs
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    filename = "/mnt/c/Users/XUKE/Desktop/hex.dump"
    f = open(filename, "w")
    close(f)
    for n = 1:100
        xyzs = Geometry.HexPoints(n, Vector{Int32}([3,3,3]), Vector{Int64}([1,1,-1]))
        Geometry.Dump(xyzs, filename, n)
    end
end



