mutable struct Point 
    index::UInt32
    type::UInt8
    coord::Vector{Int64}
    neighbors::Vector{Point}
    function Point(type::UInt8, coord::Vector{Int64})
        @assert(type==1 || type==2, "Bad point type")
        neighbors = Point[]
        new(UInt32(0), type, coord, neighbors)
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
    defects::Vector{Defect}
    defectNum::UInt32
    nstep::Int64
    function Universe(mapsize::Vector{Int64})
        map = zeros(UInt32, mapsize[1], mapsize[2], mapsize[3]) # N x half a lattice constant
        points = Point[]
        defects = Defect[]
        new(mapsize, map, points, UInt32(0), defects, UInt32(0),0)
    end
end

function Dump(universe::Universe, filename::String, mode::String)
    

function Base.push!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0, "points overlap!")
    universe.pointNum += 1
    index = universe.pointNum 
    SetIndex(universe, point, index)
    push!(universe.points, point)
end


function SetIndex(universe::Univeres, point::Point, index::UInt32)
    point.index = UInt32
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = index
end


function Base.delete!(universe::Universe, point::Point)
    universe.pointNum -= 1
    deleteat!(universe.points, point.index)
    for point in universe.points[point.index+1, end]
        SetIndex(universe, point, point.index - 1)
    end
end


function Displace!(universe::Universe, point::Point, newCoord::vector{Int64})
    universe.map[point.coord[1], point.coord[2], point[3]] = UInt32(0)
    point.coord = newCoord
    universe.map[newCoord[1], newCoord[2], newCoord[3]] = point.index
end

universe = Universe(10)

