include("geometry.jl")
using .Geometry

const DIRECTION_DIR = ([[1,1,1],
                        [1,1,-1],
                        [1,-1,1],
                        [1,-1,-1],
                        [0,0,0]])
const DEFECT_TYPE_NAME = ["SIA", "Vac"]

mutable struct Defect
    id::UInt64    # name 
    type::UInt8
    directionIndex::UInt8
    pointIndexes::Vector{UInt64}
end


mutable struct Point
    index::UInt64   # name and index (index nerver changed)
    coord::Vector{Int16}
    defect::Defect
    neighbors::Vector{Point}
    function Point(coord)
        defect = Defect(UInt64(0), UInt8(0), UInt8(0))
        neighbors = Vector{Point}[]
        new(0, coord, defect, neighbors)
    end
end

function Base.display(point::Point)
    println("$(point.index) ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
end


mutable struct Universe
    points::Vector{Point}
    defects::Vector{Defect}
    maxDefectId::UInt64
    map::Array{UInt32,3}
    mapSize::Vector{UInt16}
    function Universe(mapSize::Vector{Int64})
        mapSize = Vector{UInt16}(mapSize)
        points = Point[]
        defects = Defect[]
        map = zeros(UInt32, mapSize[1], mapSize[2], mapSize[3])
        new(points, defects, 0, map, mapSize)
    end
end

function Base.display(universe::Universe, defect::Defect)
    direction = DIRECTION_DIR[defect.directionIndex]
    println("id: $(defect.id)  type: $(DEFECT_TYPE_NAME[defect.type])  \
            direction: ($(direction[1]) $(direction[2]) $(direction[3]))")
    for index in defect.pointIndexes
        display(universe.points[index])
    end
end


function MapCoord(universe::Universe, coord::Vector{Int16})
    for i in 1:3
        if coord[i] > universe.mapSize[i]
            coord[i] -= universe.mapSize[i] 
        elseif coord[i] < 1
            coord[i] += universe.mapSize[i]
        end
    end
    coord
end

function Base.push!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    for point in points
        if universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0
            BasicInPush!(universe, point)
            MapInPush!(universe, point)
            NeighborInPush!(universe, point)
        else
            error("to do: push potin on existing point")
        end
    end
    DefectInPush!(universe, points, type, directionIndex)
    ReactInPushAndDisplace!(universe, points)
end

function ReactInPushAndDisplace!(universe::Universe, points::Vector{Point})
    defectsToMerge = Defect[]
    for point in points
        for neighbor in point.neighbors
            if neighbor.defect.type != point.defect.type
                delete!(universe, point)
                delete!(universe, neighbor)
            elseif neighbor.defect !== point.defect && neighbor.defect.type == 1 && !(neighbor.defect in defectToMerge)
                push!(defectsToMerge, neighbor.defect)
            end
        end
    end
    Merge(universe, defectsToMerge)
end

function Merge(universe::Universe, defects::Vector{Defect}) 
    
end


function DefectInPush!(universe, points, type, directionIndex)
    id = universe.maxDefectId + 1
    universe.maxDefectId = id
    pointIndexes = [p.index for p in points]
    defect = Defect(id, type, directionIndex, pointIndexes)
    [p.defect = defect for p in points]
    push!(universe.defects, defect)
end


function BasicInPush!(universe::Universe, point::Point)  # "_" means containing
    push!(universe.points, point)
    point.index = length(universe.points)
end

function MapInPush!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = point.index
end

function NeighborInPush!(universe::Universe, point::Point)
    for x in Vector{Int16}([-1, 1])
        for y in Vector{Int16}([-1, 1])
            for z in Vector{Int16}([-1, 1])
                coord = MapCoord(universe, point.coord + [x,y,z])
                neighborIndex = universe.map[coord[1], coord[2], coord[3]]
                if neighborIndex > 0
                    neighbor = universe.points[neighborIndex]
                    push!(point.neighbors, neighbor)
                    push!(neighbor.neighbors, point)
                end
            end
        end
    end
end

function AnnihilateInPush!(universe::Universe, point::Point) error("to do: annihilate") end

function Base.push!(defect::Defect, point::Point)
    push!(defect.pointIndexes, point.index)
    point.defect = defect
end





universe = Universe([100,100,100])
point1 = Point([1,1,1])
push!(universe, point1, UInt8(2), UInt8(5))