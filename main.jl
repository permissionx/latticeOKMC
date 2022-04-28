include("geometry.jl")
using .Geometry
using StatsBase
using Random
Random.seed!(1234)


const SIA_DIRECTIONS = ([Vector{Int16}([1,1,1]), 
                         Vector{Int16}([1,1,-1]), 
                         Vector{Int16}([1,-1,1]), 
                         Vector{Int16}([1,-1,-1]),
                         Vector{Int16}([0,0,0])])


const DISPLACE_DIRECTIONS = ([Vector{Int16}([1,1,1]), 
                              Vector{Int16}([1,1,-1]), 
                              Vector{Int16}([1,-1,1]), 
                              Vector{Int16}([1,-1,-1]),
                              Vector{Int16}([-1,1,1]), 
                              Vector{Int16}([-1,1,-1]), 
                              Vector{Int16}([-1,-1,1]), 
                              Vector{Int16}([-1,-1,-1])])

const DEFECT_TYPE_NAMES = ["SIA", "Vac"]


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
    print("$(point.index) $(DEFECT_TYPE_NAMES[point.defect.type]) $(point.defect.id) \
          ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
    println()
end

function Base.display(points::Vector{Point})
    println("$(length(points))-point array:")
    println("id type defect (x y z)")  # id is index
    for point in points
        display(point)
    end
    println()
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

function Base.display(universe::Universe, defect::Defect)
    direction = SIA_DIRECTIONS[defect.directionIndex]
    print("id: $(defect.id)  type: $(DEFECT_TYPE_NAMES[defect.type]) ")
    if defect.type == 1
        println("direction: ($(direction[1]) $(direction[2]) $(direction[3]))")
    else
        println()
    end
    display(universe.points[defect.pointIndexes])
end

function Base.display(universe::Universe)
    for defect in universe.defects
        display(universe, defect)
    end
end

function AlivePoints(universe::Universe)
    alivePoints = Point[]
    for defect in universe.defects
        for index in defect.pointIndexes
            push!(alivePoints, universe.points[index])
        end
    end
    return alivePoints
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
    write(file, "ITEM: ATOMS id type defect x y z\n")
    for defect in universe.defects
        for index in defect.pointIndexes
            point = universe.points[index]
            write(file, 
            "$(point.index) $(defect.type) $(defect.id) $(point.coord[1]) $(point.coord[2]) $(point.coord[3])\n")
        end
    end
    close(file)
end




function CoordInPBC(universe::Universe, coord::Vector{Int16})
    for i in 1:3
        if coord[i] > universe.mapSize[i]
            coord[i] -= universe.mapSize[i] 
        elseif coord[i] < 1
            coord[i] += universe.mapSize[i]
        end
    end
end

function CoordInPBC(universe::Universe, coords::Matrix{Int16})
    for i in 1:size(coords)[1]
        for j in 1:3
            if coords[i,j] > universe.mapSize[j]
                coords[i,j] -= universe.mapSize[j] 
            elseif coords[i,j] < 1
                coords[i,j] += universe.mapSize[j]
            end
        end
    end
end

function CoordInBCC(coord::Vector{Float64})
    # To do: should shift to the nearest lattice point, but does not matter seriously
    newCoord = Vector{Int16}(undef, 3)
    newCoord[3] = Int16(round(coord[3]))
    if iseven(newCoord[3])
        for i in 1:2
            x = floor(coord[i])
            if iseven(x)
                newCoord[i] = Int16(x)
            else
                newCoord[i] = Int16(x+1)
            end
        end
    else
        for i in 1:2
            x = floor(coord[i])
            if isodd(x)
                newCoord[i] = Int16(x)
            else
                newCoord[i] = Int16(x+1)
            end
        end
    end
    newCoord
end


function Base.push!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    # push defect by defect only 
    # vac is alwyas pushed solely
    alivePoints = Point[]
    for point in points
        isDeleted = OccupyInPushAndDisplace!(universe, point ,type)
        if !isDeleted
            BasicInPush!(universe, point)
            MapInPushAndDisplace!(universe, point)
            NeighborInPush!(universe, point)
            push!(alivePoints, point)
        end
    end
    points = alivePoints
    if length(points) > 0
        DefectInPush!(universe, points, type, directionIndex)
        ReactInPushAndDisplace!(universe, points)
    end
end



function OccupyInPushAndDisplace!(universe::Universe, point::Point, type::UInt8)
    positionIndex = universe.map[point.coord[1], point.coord[2], point.coord[3]]
    if positionIndex != 0
        # three condition:
        # vac to vac 
        # vac to sia or sia to vac
        # sia to sia
        occupiedPoint = universe.points[positionIndex]
        if occupiedPoint.defect.type != type # for vac to sia
            delete!(universe, occupiedPoint)
            return true
            # If SIA ocuupied by SIA/SIA cluster, it must be a neighbor, and thus it is goning to be merged. 
        else  # for vac to vac & SIA to SIA
            x0 = sample(Int16[1, -1])
            y0 = sample(Int16[1, -1])
            z0 = sample(Int16[1, -1])
            count = 0
            while true
                count += 1
                x = sample([x0*Int16(2), Int16(0)], Weights([1,1]))
                y = sample([y0*Int16(2), Int16(0)], Weights([1,1]))
                if x == 0 && y == 0
                    z = z0*Int16(2)
                else
                    z = sample([z0*Int16(2), Int16(0)], Weights([1,1]))
                end
                if x == 2 && y == 2 && z == 2
                    x = Int16(1)
                    y = Int16(1)
                    z = Int16(1)
                end
                point.coord = [point.coord[1]+x, point.coord[2]+y, point.coord[3]+z]
                CoordInPBC(universe, point.coord)
                if universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0
                    break
                end
                @assert(count < 1000, "Too many iterations in OccupyInPushAndDisplace!")
            end
        end
    end
    return false
end

function BasicInPush!(universe::Universe, point::Point)  
    push!(universe.points, point)
    point.index = length(universe.points)
    universe.pointNum += 1
end

function MapInPushAndDisplace!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = point.index
end

function NeighborInPush!(universe::Universe, point::Point)
    for x in Vector{Int16}([-1, 1])
        for y in Vector{Int16}([-1, 1])
            for z in Vector{Int16}([-1, 1])
                coord = point.coord + [x,y,z]
                CoordInPBC(universe, coord)
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

function DefectInPush!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    id = universe.maxDefectId + 1
    universe.maxDefectId = id
    pointIndexes = [point.index for point in points]
    defect = Defect(id, type, directionIndex, pointIndexes)
    [point.defect = defect for point in points]
    push!(universe.defects, defect)
end


function ReactInPushAndDisplace!(universe::Universe, points::Vector{Point})
    

    defect = points[1].defect
    defectsToMerge = Defect[defect]
    for point in points
        deleteNeighbors = Point[]
        for neighbor in point.neighbors
            if neighbor.defect.type != defect.type
                push!(deleteNeighbors, neighbor)
            elseif neighbor.defect.type == 1 && neighbor.defect !== defect &&  !(neighbor.defect in defectsToMerge)
                push!(defectsToMerge, neighbor.defect)
            end
        end
        if length(deleteNeighbors) > 0
            neighbor = sample(deleteNeighbors)
            delete!(universe, neighbor)
            delete!(universe, point)
        end
    end

    

    if length(defectsToMerge) > 1 
        Merge!(universe, defectsToMerge)
    end
end


function Merge!(universe::Universe, defects::Vector{Defect}) 
    # merge defects to defects[1]
    MergeBasic!(universe, defects)
    Arrange!(universe, defects[1])
end


function MergeBasic!(universe::Universe, defects::Vector{Defect})
    defect = defects[1]
    for i in 2:length(defects)
        newPointIndexes = defects[i].pointIndexes
        defect.pointIndexes = vcat(defect.pointIndexes, newPointIndexes)
        for index in newPointIndexes
            universe.points[index].defect = defect
        end
        deleteat!(universe.defects, findfirst(x -> x === defects[i], universe.defects))   
    end
end

function Arrange!(universe::Universe, defect::Defect)
    aveCoord = Vector{Int16}([0,0,0])
    for pointIndex in defect.pointIndexes
        point = universe.points[pointIndex]
        aveCoord += point.coord
    end
    aveCoord = CoordInBCC(aveCoord / length(defect.pointIndexes))
    coords = HexPoints(length(defect.pointIndexes), aveCoord, SIA_DIRECTIONS[defect.directionIndex])
    points = [universe.points[pointIndex] for pointIndex in defect.pointIndexes]
    displace!(universe, points, coords)
end

function Base.delete!(universe::Universe, point::Point)
    BasicInDelete!(universe, point)
    MapInDelete!(universe, point)
    NeighborInDelete!(universe, point)
    DefectInDelete!(universe, point)
end

function BasicInDelete!(universe::Universe, point::Point) 
    universe.pointNum -= 1
end

function MapInDelete!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
end

function NeighborInDelete!(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(x -> x === point, neighbor.neighbors))
    end
end

function DefectInDelete!(universe::Universe, point::Point)
    defect = point.defect
    if length(defect.pointIndexes) == 1
        deleteat!(universe.defects, findfirst(x -> x === defect, universe.defects))
    else
        deleteat!(defect.pointIndexes, findfirst(x -> x == point.index, defect.pointIndexes))
    end
end


function displace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int16}) 
    CleanMapInDisplace!(universe, points)
    CoordInPBC(universe, newCoords)
    alivePoints = Point[]
    for point in points
        BasicInDisplace!(points, newCoords)
        isDeleted = OccupyInPushAndDisplace!(universe, point, point.defect.type)
        if !isDeleted
            MapInPushAndDisplace!(universe, point)
            NeighborInDisplace!(universe, point)
            push!(alivePoints, point)
        end
    end
    points = alivePoints
    ReactInPushAndDisplace!(universe, points)
end

function CleanMapInDisplace!(universe::Universe, points::Vector{Point})
    for point in points
        universe.map[point.coord[1], point.coord[2], point.coord[3]] = UInt32(0)
    end
end

function BasicInDisplace!(points::Vector{Point}, newCoords::Matrix{Int16})
    for i in 1:length(points)
        points[i].coord = newCoords[i, :]
    end
end


function NeighborInDisplace!(universe::Universe, point::Point)
    point.neighbors = Point[]
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(x -> x === point, neighbor.neighbors))
    end
    for x in Vector{Int16}([-1, 1])
        for y in Vector{Int16}([-1, 1])
            for z in Vector{Int16}([-1, 1])
                coord = point.coord + [x,y,z]
                CoordInPBC(universe, coord)
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


mapSize = Vector{UInt16}([300,300,300])
universe = Universe(mapSize)
fileName = "/mnt/c/Users/xuke/Desktop/test.dump"
RefreshFile(fileName)

function test1!(universe::Universe)
    point1 = Point(Vector{Int16}([11,11,11]))
    point2 = Point(Vector{Int16}([12,12,12]))
    point3 = Point(Vector{Int16}([13,13,13]))
    point4 = Point(Vector{Int16}([14,14,14]))
    point5 = Point(Vector{Int16}([15,15,15]))
    #Base.push!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    push!(universe, [point1, point2, point3, point4, point5], UInt8(2), UInt8(1))
    #universe.nStep += 1
    #Dump(universe, fileName, "a")
    #universe.nStep += 1
    #Dump(universe, fileName, "a")
    # To do: add point on exsiting point

    point6 = Point(Vector{Int16}([13,13,13]))
    push!(universe, [point6], UInt8(2), UInt8(1))
    Dump(universe, fileName, "a")
end

function test2!(universe::Universe)
    for i in 1:1000
        universe.nStep += 1
        point = Point(Vector{Int16}([150,150,150]))
        push!(universe, [point], UInt8(1), UInt8(1))
        if universe.nStep % 10 == 0
            println("step: ", universe.nStep)
            Dump(universe, fileName, "a")
        end
    end
end
test2!(universe)

