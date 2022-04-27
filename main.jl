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
        pointIndexes = UInt64[]
        defect = Defect(UInt64(0), UInt8(0), UInt8(0), pointIndexes)
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
    # push defect by defect only 
    for point in points
        if universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0
            BasicInPush!(universe, point)
            MapInPush!(universe, point)
            NeighborInPush!(universe, point)
        else
            error("to do: push potin on existing point")
        end
    end
    defect = DefectInPush!(universe, points, type, directionIndex)
    ReactInPushAndDisplace!(universe, points, defect)
end


function ReactInPushAndDisplace!(universe::Universe, points::Vector{Point}, defect::Defect)
    defectsToMerge = Defect[defect]
    for point in points
        for neighbor in point.neighbors
            if neighbor.defect.type != defect.type
                delete!(universe, point)
                delete!(universe, neighbor)
            elseif neighbor.defect.type == 1 && neighbor.defect !== defect &&  !(neighbor.defect in defectToMerge)
                push!(defectsToMerge, neighbor.defect)
            end
        end
    end
    if length(defectsToMerge) > 1 
        Merge!(universe, defectsToMerge)
    end
end


function Merge!(universe::Universe, defects::Vector{Defect}) 
    # merge defects to defects[1]
    defect = MergeBasic!(universe, defects)
    Rearrange!(universe, defect)
end


function MergeBasic!(universe::Universe, defects::Vector{Defect})
    for i in 2:length(defects)
        newPointIndexes = defects[i].pointIndexes
        defects[1].pointIndexes = vcat(defects[1].pointIndexes, newPointIndexes)
        for index in newPointIndexes
            universe.points[index].defect = defects[1]
        end
        deleteat!(universe.defects, findfirst(d -> d === defect, universe.defects))    
    end
end



function Base.delete!(universe::Universe, point::Point)
    BasicInDelete!(universe, point)
    MapInDelete!(universe, point)
    NeighborInDelete!(universe, point)
    DefectInDelete!(universe, point)
end

function BasicDelete!(universe::Universe, point::Point) end

function MapInDelete!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
end

function NeighborInDelete!(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(p -> p === point, neighbor.neighbors))
    end
end

function DefectInDelete!(universe::Universe, point::Point)
    defect = point.defect
    if length(defect.pointIndexes) == 1
        deleteat!(universe.defects, findfirst(d -> d === defect, universe.defects))
    else
        deleteat!(defect.pointIndexes, findfirst(pIndex -> pIndex == point.index, defect.pointIndexes))
    end
end

function BasicInDelete!(universe::Universe, point::Point) end



function Rearrange!(universe::Universe, defect::Defect)
    aveCoord = Vector{Int16}([0,0,0])
    for pIndex in defect.directionIndex
        p = universe.points[pIndex]
        aveCoord += p.coord
    end
    aveCoord = Vector{Int16}(round.(aveCoord / length(defect.directionIndex)))
    coords = HexPoints(length(defect.directionIndex), aveCoord, DIRECTION_DIR[defect.directionIndex])
    points = [universe.points[pIndex] for pIndex in defect.pointIndexes]
    displace!(universe, points, coords)
end

function displace!(universe::Universe, points::Vector{Point}, newCoords::Vector{Int16}) 
    BasicInDisplace!(universe, points, newCoord)
    MapInDisplace!(universe, points, newCoord)
    NeighborInDisplace!(universe, points, newCoord)
    DefectInDisplace!(universe, points, newCoord)
end

function BasicInDisplace!(universe::Universe, points::Vector{Point}, newCoord::Vector{Int16})
    for p in points
        p.coord = newCoord
    end
end

function MapInDisplace!(universe::Universe, points::Vector{Point}, newCoord::Vector{Int16})
    for p in points
        universe.map[p.coord[1], p.coord[2], p.coord[3]] = 0
    end
    for p in points
        universe.map[newCoord[1], newCoord[2], newCoord[3]] = p.index
    end
end

function NeighborInDisplace!(universe::Universe, point::Points, newCoord::Vector{Int64}) 
    
end




function DefectInPush!(universe, points, type, directionIndex)
    id = universe.maxDefectId + 1
    universe.maxDefectId = id
    pointIndexes = [p.index for p in points]
    defect = Defect(id, type, directionIndex, pointIndexes)
    [p.defect = defect for p in points]
    push!(universe.defects, defect)
    defect
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



universe = Universe([100,100,100])
point1 = Point([1,1,1])
push!(universe, point1, UInt8(2), UInt8(5))

