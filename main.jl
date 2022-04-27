include("geometry.jl")
using .Geometry


const DIRECTION_DIR = ([[1,1,1],
                        [1,1,-1],
                        [1,-1,1],
                        [1,-1,-1],
                        [0,0,0]])
const DEFECT_TYPE_NAME = ["SIA", "Vac"]


mutable struct Defect
    id::UInt32    # name 
    type::UInt8
    directionIndex::UInt8
    pointIndexes::Vector{UInt32}
end


mutable struct Point
    index::UInt32   # name and index (index nerver changed)
    coord::Vector{Int16}
    defect::Defect
    neighbors::Vector{Point}
    function Point(coord)
        pointIndexes = UInt32[]
        defect = Defect(UInt32(0), UInt8(0), UInt8(0), pointIndexes)
        neighbors = Vector{Point}[]
        new(0, coord, defect, neighbors)
    end
end


function Base.display(point::Point)
    println("$(point.index) ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
end


mutable struct Universe
    points::Vector{Point}
    pointNum::UInt32
    defects::Vector{Defect}
    maxDefectId::UInt32
    map::Array{UInt32,3}
    mapSize::Vector{UInt16}
    nStep::Int64
    function Universe(mapSize::Vector{UInt16})
        mapSize = Vector{UInt16}(mapSize)
        points = Point[]
        defects = Defect[]
        map = zeros(UInt32, mapSize[1], mapSize[2], mapSize[3])
        new(points, UInt32(0), defects, UInt32(0), map, mapSize, 0)
    end
end

function RefreshFile(fileName::String)
    file = open(fileName, "w")
    close(file)
end

function Dump(universe::Universe, fileName::String, mode::String)
    file = open(fileName, mode)
    write(file, "ITEM: TIMESTEP\n")
    write(file, "$(universe.nStep)\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
    write(file, "$(universe.pointNum)\n")
    write(file, "ITEM: BOX BOUNDS pp pp pp\n")
    write(file, "1 $(universe.mapSize[1]+1)\n")
    write(file, "1 $(universe.mapSize[2]+1)\n")
    write(file, "1 $(universe.mapSize[3]+1)\n")
    write(file, "ITEM: ATOMS id type x y z\n")
    for i = 1:length(universe.defects)
        for j = 1:length(defect.pointIndexes)
            point = universe.points[j]
            write(file, 
            "$(point.index) $(defect.type) $(point.coord[1]) $(point.coord[2]) $(point.coord[3])\n")
        end
    end
    close(file)
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
    ReactInPushAndDisplace!(universe, points)
end

function BasicInPush!(universe::Universe, point::Point)  # "_" means containing
    push!(universe.points, point)
    point.index = length(universe.points)
    universe.pointNum += 1
end

function DefectInPush!(universe, points, type, directionIndex)
    id = universe.maxDefectId + 1
    universe.maxDefectId = id
    pointIndexes = [point.index for point in points]
    defect = Defect(id, type, directionIndex, pointIndexes)
    [point.defect = defect for point in points]
    push!(universe.defects, defect)
    defect
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

function ReactInPushAndDisplace!(universe::Universe, points::Vector{Point})
    defect = points[1].defect
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
    Arrange!(universe, defect)
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

function Arrange!(universe::Universe, defect::Defect)
    aveCoord = Vector{Int16}([0,0,0])
    for pIndex in defect.directionIndex
        point = universe.points[pIndex]
        aveCoord += point.coord
    end
    aveCoord = Vector{Int16}(round.(aveCoord / length(defect.directionIndex)))
    coords = HexPoints(length(defect.directionIndex), aveCoord, DIRECTION_DIR[defect.directionIndex])
    points = [universe.points[pIndex] for pIndex in defect.pointIndexes]
    displace!(universe, points, coords)
end

function Base.delete!(universe::Universe, point::Point)
    BasicInDelete!(universe, point)
    MapInDelete!(universe, point)
    NeighborInDelete!(universe, point)
    DefectInDelete!(universe, point)
end

function BasicDelete!(universe::Universe, point::Point) 
    universe.pointNum -= 1
end

function MapInDelete!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
end

function NeighborInDelete!(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(point -> point === point, neighbor.neighbors))
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



function displace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int16}) 
    newCoords = [MapCoord(universe, newCoords[i,:]) for i in size(newCoords)[1]]
    BasicAndMapInDisplace!(universe, points, newCoords)
    NeighborInDisplace!(universe, points)
    ReactInPushAndDisplace!(universe, points)
end



function BasicAndMapInDisplace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int16})
    for i in 1:length(points)
        point = points[i]
        universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
        point.coord = newCoords[i,:]
        universe.map[point.coord[1], point.coord[2], point.coord[3]] = point.index
    end
end

function NeighborInDisplace!(universe::Universe, points::Vector{Point}) 
    for point in points
        oldNeighborIndexes = [point.index for point in point.neighbors]
        neighborIndexes = UInt32[]
        for x in Vector{Int16}([-1, 1])
            for y in Vector{Int16}([-1, 1])
                for z in Vector{Int16}([-1, 1])
                    coord = MapCoord(universe, point.coord + [x,y,z])
                    neighborIndex = universe.map[coord[1], coord[2], coord[3]]
                    if neighborIndex > 0
                        push!(neighborIndexes, neighborIndex)
                    end
                end
            end
        end
        if !issetequal(oldNeighborIndexes, neighborIndexes)
            dropedNeighborIndexes = setdiff(neighborIndexesBefore, neighborIndexes)
            newNeighborIndexes = setdiff(neighborIndexes, neighborIndexesBefore)
            for dropedNeighborIndex in dropedNeighborIndexes
                dropedNeighbor = universe.points[dropedNeighborIndex]
                deleteat!(dropedNeighbor.neighbors, 
                          findfirst(point -> point === point, dropedNeighbor.neighbors))
            end
            for newNeighborIndex in newNeighborIndexes
                newNeighbor = universe.points[newNeighborIndex]
                push!(newNeighbor.neighbors, point)
            end
            point.neighbors = [universe.points[index] for index in neighborIndexes]
        end
    end
end



mapSize = Vector{UInt16}([10,10,10])
universe = Universe(mapSize)
fileName = "/mnt/c/Users/buaax/Desktop/test.dump"
RefreshFile(fileName)
Dump(universe, fileName, "a")
