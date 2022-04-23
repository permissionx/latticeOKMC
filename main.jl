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
    function Universe(mapsize::Vector{Int64})
        map = zeros(UInt32, mapsize[1], mapsize[2], mapsize[3]) # N x half a lattice constant
        points = Point[]
        defects = Defect[]
        new(mapsize, map, points, UInt32(0), defects, UInt32(0))
    end
end


function FeindNeighbors!(universe, point)
    neighbors = Point[]
    for x in [-1,1]
        for y in [-1,1]
            for z in [-1,1]
                neighbor = universe.points[universe.map[x,y,z]]
                push!(neighbors, neighbor)
            end
        end
    end
end


function Base.push!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0, "points overlap!")
    map.content[point.coord[1], point.coord[2], point.coord[3]] = point.index
    universe.pointNum += 1
    index = universe.pointNum 
    point.index = index
    push!(universe.points, point)
    FindNeighbor!(universe, point)
end


function Base.delete!(map::Map, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]])


function Base.delete!(universe::Universe, point::Point)
    universe.pointNum -= 1
    deleteat!(universe.points, point.index)
    for point in universe.points[point.index+1, end]
        point.index -= 1
    end
    delete!(univeres.map, point)
end

universe = Universe(10)
