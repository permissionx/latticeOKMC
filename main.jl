include("geometry.jl")
using .Geometry
using StatsBase
using Random
using Distributions
Random.seed!(1234)


const SIA_DIRECTIONS = ([Vector{Int32}([1,1,1]), 
                         Vector{Int32}([1,1,-1]), 
                         Vector{Int32}([1,-1,1]), 
                         Vector{Int32}([1,-1,-1])])


const DISPLACE_DIRECTIONS = ([Vector{Int32}([1,1,1]), 
                              Vector{Int32}([1,1,-1]), 
                              Vector{Int32}([1,-1,1]), 
                              Vector{Int32}([1,-1,-1]),
                              Vector{Int32}([-1,1,1]), 
                              Vector{Int32}([-1,1,-1]), 
                              Vector{Int32}([-1,-1,1]), 
                              Vector{Int32}([-1,-1,-1])])

                              
const DEFECT_TYPE_NAMES = ["SIA", "Vac"]


mutable struct Defect
    id::UInt32    # name 
    type::UInt8
    directionIndex::UInt8
    pointIndexes::Vector{UInt32}
end


mutable struct Point
    index::UInt32   # name and index (index nerver changed)
    coord::Vector{Int32}
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
    mapSize::Vector{UInt32}
    nStep::Int64
    function Universe(mapSize::Vector{UInt32})
        mapSize = Vector{UInt32}(mapSize)
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


function PBCCoord!(universe::Universe, coord::Vector{Int32})
    for i in 1:3
        if coord[i] > universe.mapSize[i]
            coord[i] -= universe.mapSize[i] 
        elseif coord[i] < 1
            coord[i] += universe.mapSize[i]
        end
    end
end


function PBCCoord!(universe::Universe, coords::Matrix{Int32})
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
            x0 = sample(Int32[1, -1])
            y0 = sample(Int32[1, -1])
            z0 = sample(Int32[1, -1])
            count = 0
            while true
                count += 1
                x = sample([x0*Int32(2), Int32(0)], Weights([1,1]))
                y = sample([y0*Int32(2), Int32(0)], Weights([1,1]))
                if x == 0 && y == 0
                    z = z0*Int32(2)
                else
                    z = sample([z0*Int32(2), Int32(0)], Weights([1,1]))
                end
                if x == 2 && y == 2 && z == 2
                    x = Int32(1)
                    y = Int32(1)
                    z = Int32(1)
                end
                point.coord = [point.coord[1]+x, point.coord[2]+y, point.coord[3]+z]
                PBCCoord!(universe, point.coord)
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
    for x in Vector{Int32}([-1, 1])
        for y in Vector{Int32}([-1, 1])
            for z in Vector{Int32}([-1, 1])
                coord = point.coord + [x,y,z]
                PBCCoord!(universe, coord)
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
        isDeleted = false
        for neighbor in point.neighbors
            if neighbor.defect.type != defect.type
                push!(deleteNeighbors, neighbor)
                isDeleted = true
            elseif neighbor.defect.type == 1 && !(isDeleted) && neighbor.defect !== defect &&  !(neighbor.defect in defectsToMerge)
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
    defect = MergeBasic!(universe, defects)
    Arrange!(universe, defect)
end


function MergeBasic!(universe::Universe, defects::Vector{Defect})
    _,index = findmax(x->length(x.pointIndexes), defects)
    largeDefect = defects[index]
    for defect in defects
        if !(defect === largeDefect)
            newPointIndexes = defect.pointIndexes
            largeDefect.pointIndexes = vcat(largeDefect.pointIndexes, newPointIndexes)
            for index in newPointIndexes
                universe.points[index].defect = largeDefect
            end
            deleteat!(universe.defects, findfirst(x -> x === defect, universe.defects))   
        end
    end
    largeDefect
end


function Arrange!(universe::Universe, defect::Defect)
    aveCoord = Vector{Int64}([0,0,0])
    for pointIndex in defect.pointIndexes
        point = universe.points[pointIndex]
        aveCoord += point.coord
    end
    aveCoord = aveCoord / length(defect.pointIndexes)
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


function displace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int32}) 
    CleanMapInDisplace!(universe, points)
    PBCCoord!(universe, newCoords)
    alivePoints = Point[]  
    BasicInDisplace!(points, newCoords)  
    for point in points
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


function BasicInDisplace!(points::Vector{Point}, newCoords::Matrix{Int32})
    for i in 1:length(points)
        points[i].coord = newCoords[i, :]
    end
end


function NeighborInDisplace!(universe::Universe, point::Point)
    point.neighbors = Point[]
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(x -> x === point, neighbor.neighbors))
    end
    for x in Vector{Int32}([-1, 1])
        for y in Vector{Int32}([-1, 1])
            for z in Vector{Int32}([-1, 1])
                coord = point.coord + [x,y,z]
                PBCCoord!(universe, coord)
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


function test1()
    mapSize = Vector{UInt32}([300,300,300])
    universe = Universe(mapSize)
    fileName = "/mnt/c/Users/xuke/Desktop/test1.dump"
    RefreshFile(fileName)
    point1 = Point(Vector{Int32}([11,11,11]))
    point2 = Point(Vector{Int32}([12,12,12]))
    point3 = Point(Vector{Int32}([13,13,13]))
    point4 = Point(Vector{Int32}([14,14,14]))
    point5 = Point(Vector{Int32}([15,15,15]))
    #Base.push!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    push!(universe, [point1, point2, point3, point4, point5], UInt8(2), UInt8(1))
    universe.nStep += 1
    Dump(universe, fileName, "a")
    universe.nStep += 1
    Dump(universe, fileName, "a")    
    point6 = Point(Vector{Int32}([13,13,13]))
    push!(universe, [point6], UInt8(2), UInt8(1))
    Dump(universe, fileName, "a")
    universe
end


function test2()
    mapSize = Vector{UInt32}([300,300,300])
    universe = Universe(mapSize)
    fileName = "/mnt/c/Users/xuke/Desktop/test2.dump"
    RefreshFile(fileName)
    for i in 1:1000
        universe.nStep += 1
        println("step: ", universe.nStep)
        point = Point(Vector{Int32}([150,150,150]))
        push!(universe, [point], UInt8(1), UInt8(1))
        if universe.nStep % 10 == 0
            Dump(universe, fileName, "a")
        end
    end
    universe
end


function test3()
    mapSize = Vector{UInt32}([300,300,300])
    universe = Universe(mapSize)
    fileName = "/mnt/c/Users/xuke/Desktop/test3.dump"
    RefreshFile(fileName)
    for i in 1:10000
        universe.nStep += 1
        println("step: ", universe.nStep)
        coord = rand(Normal(150, 15), 3)
        coord = Geometry.CoordInBCC(coord)
        PBCCoord!(universe, coord)
        point = Point(coord)
        type = sample(UInt8(1):UInt8(2))
        direction = sample(UInt8(1):UInt8(4))
        push!(universe, [point], type, direction)
        if universe.nStep % 100 == 0
            Dump(universe, fileName, "a")
        end
    end
    universe
end

universe = test3()

