mutable struct Point 
    id::Int64
    index::UInt32
    type::UInt8
    coord::Vector{Int64}
    neighbors::Vector{Point}
    function Point(type::UInt8, coord::Vector{Int64})
        @assert(type==1 || type==2, "Bad point type")
        neighbors = Point[]
        new(0, UInt32(0), type, coord, neighbors)
    end
end

const pointTypeNames = ["SIA", "Vac"]
function Base.display(point::Point)
    print("$(point.id) $(point.index) $(pointTypeNames[point.type]) \
          ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
end

function Base.display(points::Vector{Point})
    println("$(length(points))-point universe:")
    println("id index type (x y z)")
    for point in points
        display(point)
        println()
    end
end



mutable struct Defect
    index::UInt32
    type::UInt8
    points::Vector{Point}
    direction::Vector{Int64}
    function Defect()
        points = Point[]
        new(UInt32(0), UInt8(0), points, [0,0,0])
    end
end
 
mutable struct Universe
    mapsize::Vector{Int64}
    map::Array{UInt32, 3}
    points::Vector{Point}
    pointNum::UInt32
    maxID::Int64
    defects::Vector{Defect}
    defectNum::UInt32
    nstep::Int64
    function Universe(mapsize::Vector{Int64})
        map = zeros(UInt32, mapsize[1], mapsize[2], mapsize[3]) # N x half a lattice constant
        points = Point[]
        defects = Defect[]
        new(mapsize, map, points, UInt32(0), 0, defects, UInt32(0),0)
    end
end

# dump the map in lammps dump format
function Dump(universe::Universe, filename::String, mode::String)
    file = open(filename, mode)
    write(file, "ITEM: TIMESTEP\n")
    write(file, "$(universe.nstep)\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
    write(file, "$(universe.pointNum)\n")
    write(file, "ITEM: BOX BOUNDS pp pp pp\n")
    write(file, "1 $(universe.mapsize[1])\n")
    write(file, "1 $(universe.mapsize[2])\n")
    write(file, "1 $(universe.mapsize[3])\n")
    write(file, "ITEM: ATOMS id type x y z\n")
    for i=1:universe.pointNum
        point = universe.points[i]
        write(file, "$(point.id) $(point.type) $(point.coord[1]) $(point.coord[2]) $(point.coord[3])\n")
    end
    close(file)
end

using Crayons
function Base.display(universe::Universe, zorder::Int64)
    map2D = universe.map[:,:,zorder]
    for i=1:universe.mapsize[1]
        for j=1:universe.mapsize[2]
            if isodd(zorder)
                if isodd(i) && isodd(j)
                    if map2D[i,j] > 0
                        color = :light_red
                    else
                        color = :default
                    end
                    print(Crayon(foreground=color, bold=true), map2D[i,j], " ")
                else
                    print(Crayon(foreground=(105,105,105), bold=true),"X ")
                end
            else 
                if iseven(i) && iseven(j)
                    if map2D[i,j] > 0
                        color = :light_red
                    else
                        color = :default
                    end
                    print(Crayon(foreground=color, bold=true), map2D[i,j], " ")
                else
                    print(Crayon(foreground=(105,105,105), bold=true),"X ")
                end
            end
        end
        println()
    end
end


using Printf
function EmojiDisplay(universe::Universe, zorder::Int64)
    map2D = universe.map[:,:,zorder]
    for i=1:universe.mapsize[1]
        for j=1:universe.mapsize[2]
            if isodd(zorder)
                if isodd(i) && isodd(j)
                    if map2D[i,j] == 0
                        @printf "%-1.1s" "⭕️"
                    else
                        @printf "%-1.1s" "🚀"
                    end
                else
                    @printf "%-1.1s" "❌"
                end
            else 
                if iseven(i) && iseven(j)
                    if map2D[i,j] == 0
                        @printf "%-1.1s" "⭕️"
                    else
                        @printf "%-1.1s" "🚀"
                    end
                else
                    @printf "%-1.1s" "❌"
                end
            end
        end
        println()
    end
end



function RefreshFile(filename::String)
    file = open(filename, "w")
    close(file)
end


function Base.push!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0, "points overlap!")
    PushBasic!(universe, point)
    PushMap!(universe, point)
    PushNeighbors!(universe, point)
end

function PushBasic!(universe::Universe, point::Point)
    universe.pointNum += 1
    universe.maxID += 1
    point.index = universe.pointNum
    point.id = universe.maxID
    push!(universe.points, point)
end

function PushMap!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = point.index
end

function PushNeighbors!(universe::Universe, point::Point)
    for i in [-1,1]
        for j in [-1,1]
            for k in [-1,1]
                coord = point.coord + [i,j,k]
                try
                    if universe.map[coord[1], coord[2], coord[3]] > 0
                        neighbor = universe.points[universe.map[coord[1], coord[2], coord[3]]]
                        push!(point.neighbors, neighbor)
                        push!(neighbor.neighbors, point)
                    end
                catch BoundsError
                    continue
                end
            end
        end
    end
end

function Base.delete!(universe::Universe, point::Point)
    DeleteNeighbors(universe, point)
    DeleteBasic!(universe, point)
    DeleteMap!(universe, point)
end

function DeleteBasic!(universe::Universe, point::Point)
    universe.pointNum -= 1
    deleteat!(universe.points, point.index)
    for point in universe.points[point.index:end]
        point.index -= 1
    end
end

function DeleteMap!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] > 0, "point not in map!")
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
    for p in universe.points[point.index:end]
        universe.map[p.coord[1], p.coord[2], p.coord[3]] -= 1
    end
end

function DeleteNeighbors(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(p->p.id==point.id,neighbor.neighbors))
    end
end


function Displace!(universe::Universe, point::Point, newCoord::Vector{Int64})
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] > 0, "point not in map!")
    @assert(universe.map[newCoord[1], newCoord[2], newCoord[3]] == 0, "points overlap!")
    DisplaceBasic!(universe, point, newCoord)
end

function DisplaceBasic!(universe::Universe, point::Point, newCoord::Vector{Int64})
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = UInt32(0)
    point.coord = newCoord
    universe.map[newCoord[1], newCoord[2], newCoord[3]] = point.index
end


universe = Universe([10,10,10])
point1 = Point(UInt8(1), [3,3,3])
point2 = Point(UInt8(1), [4,4,4])
point3 = Point(UInt8(1), [2,2,2])

filename = "/mnt/c/Users/buaax/Desktop/test.dump"
RefreshFile(filename)
push!(universe, point1)
push!(universe, point2)
push!(universe, point3)
Dump(universe, filename, "a")