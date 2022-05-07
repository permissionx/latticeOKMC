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
    write(file, "ITEM: ATOMS id type defect object x y z\n")
    for defect in universe.defects
        for index in defect.pointIndexes
            point = universe.points[index]
            write(file, 
            "$(point.index) $(defect.type) $(defect.index) $(point.object.index) $(point.coord[1]) $(point.coord[2]) $(point.coord[3])\n")
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
    alivePoints = Point[]
    for point in points
        isDeleted = OccupyInPushAndDisplace!(universe, point ,type)
        if !isDeleted
            BasicInPush!(universe, point)
            MapInPushAndDisplace!(universe, point)
            NeighborInPush!(universe, point)
            if type === UInt8(2)
                ObjectInPushVac(universe, point)
            end
            push!(alivePoints, point)
        end
    end
    points = alivePoints
    if length(points) > 0
        defect = DefectInPush!(universe, points, type, directionIndex)
        if type === UInt8(1)
            ObjectInPushSia(universe, points, defect)
        end
        ReactInPushAndDisplace!(universe, points)
    end
end



function PushRefreshList!(universe::Universe, object::Object)
    if !(object in universe.objectsToRefresh)
        push!(universe.objectsToRefresh, object)
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
        if occupiedPoint.type != type # for vac to sia
            delete!(universe, occupiedPoint)
            return true
            # If Sia ocuupied by Sia/Sia cluster, it must be a neighbor, and thus it is goning to be merged. 
        else  # for vac to vac & Sia to Sia
            x0 = sample(Int32[1, -1])
            y0 = sample(Int32[1, -1])
            z0 = sample(Int32[1, -1])
            displaceVectors = [Int32[2*x0,0,0],
                               Int32[0,2*y0,0],
                               Int32[0,0,2*z0],
                               Int32[1*z0,1*y0,1*z0]]
            count = 0
            while true
                count += 1
                displaceVector = sample(displaceVectors)
                point.coord += displaceVector
                PBCCoord!(universe, point.coord)
                if universe.map[point.coord[1], point.coord[2], point.coord[3]] === UInt32(0)
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
    for neighborVector in NEIGHBOR_VECTORS
        coord = point.coord + neighborVector
        PBCCoord!(universe, coord)
        neighborIndex = universe.map[coord[1], coord[2], coord[3]]
        if neighborIndex > 0
            neighbor = universe.points[neighborIndex]
            push!(point.neighbors, neighbor)
            push!(neighbor.neighbors, point)
            PushRefreshList!(universe, neighbor.object)
        end
    end
end


function DefectInPush!(universe::Universe, points::Vector{Point}, type::UInt8, directionIndex::UInt8)
    index = universe.maxDefectIndex + 1
    universe.maxDefectIndex = index
    pointIndexes = [point.index for point in points]
    defect = Defect(index, type, directionIndex, pointIndexes)
    for point in points
        point.defect = defect
        point.type = defect.type
    end
    push!(universe.defects, defect)
    defect
end


function ObjectInPushVac(universe::Universe, point::Point)
    universe.maxObjectIndex += 1
    index = universe.maxObjectIndex
    pointIndexes = [point.index]
    object = Object(index, UInt8(2), UInt8(1), pointIndexes)
    point.object = object
    push!(universe.objects, object)
    PushRefreshList!(universe, object)
end


function ObjectInPushSia(universe::Universe, points::Vector{Point}, defect::Defect)
    universe.maxObjectIndex += 1
    index = universe.maxObjectIndex
    object = Object(index, UInt8(1), defect.directionIndex, defect.pointIndexes)
    for point in points
        point.object = object
    end
    push!(universe.objects, object)
    PushRefreshList!(universe, object)
end


function ReactInPushAndDisplace!(universe::Universe, points::Vector{Point})
    defect = points[1].defect
    defectsToMerge = Defect[defect]
    for point in points
        deleteNeighbors = Point[]
        for neighbor in point.neighbors
            if neighbor.type != defect.type
                push!(deleteNeighbors, neighbor)
            end
        end
        if length(deleteNeighbors) > 0
            neighbor = sample(deleteNeighbors)
            delete!(universe, neighbor)
            delete!(universe, point)
            continue
        end
        for neighbor in point.neighbors
            if neighbor.defect !== defect 
                if !(neighbor.defect in defectsToMerge)
                    push!(defectsToMerge, neighbor.defect)
                end
            end
        end
    end
    if length(defectsToMerge) > 1
        Merge!(universe, defectsToMerge)
    end
end

function Merge!(universe::Universe, defects::Vector{Defect}) 
    # merge defects to defects[1]
    defect = DefectAndObjectInMerge!(universe, defects)
    if defect.type === UInt8(1)
        Arrange!(universe, defect)
    end
end


function DefectAndObjectInMerge!(universe::Universe, defects::Vector{Defect})
    _,index = findmax(x->length(x.pointIndexes), defects)
    largeDefect = defects[index]
    type = largeDefect.type
    if type === UInt8(1)
        largeObject = universe.points[largeDefect.pointIndexes[1]].object
    end
    for defect in defects
        if !(defect === largeDefect)
            newPointIndexes = defect.pointIndexes
            largeDefect.pointIndexes = vcat(largeDefect.pointIndexes, newPointIndexes)
            object = universe.points[defect.pointIndexes[1]].object
            for index in newPointIndexes
                point = universe.points[index]
                point.defect = largeDefect
                if type === UInt8(1)
                    point.object = largeObject
                end
            end
            deleteat!(universe.defects, findfirst(x -> x === defect, universe.defects))  
            if type === UInt8(1)
                delete!(universe, object) 
            end 
        end
    end
    if type === UInt8(1)
        largeObject.pointIndexes = largeDefect.pointIndexes
        PushRefreshList!(universe, largeObject)
    end
    largeDefect
end

function Base.delete!(universe::Universe, object::Object)
    deleteat!(universe.objects, findfirst(x -> x === object, universe.objects))
    if object in universe.objectsToRefresh
        deleteat!(universe.objectsToRefresh, findfirst(x -> x === object, universe.objectsToRefresh))
    end
end


function Arrange!(universe::Universe, defect::Defect)
    aveCoord = Vector{Int64}([0,0,0])
    coord1 = universe.points[defect.pointIndexes[1]].coord
    for pointIndex in defect.pointIndexes
        point = universe.points[pointIndex]
        coord = copy(point.coord)
        for i in 1:3
            if coord[i] - coord1[i] > mapSize[i] / 2
                coord[i] -= mapSize[i]
            elseif coord[i] - coord1[i] < -mapSize[i] / 2
                coord[i] += mapSize[i]
            end
        end
        aveCoord += coord
    end
    aveCoord = aveCoord / length(defect.pointIndexes)
    coords = HexPoints(length(defect.pointIndexes), aveCoord, Sia_DIRECTIONS[defect.directionIndex])
    points = [universe.points[pointIndex] for pointIndex in defect.pointIndexes]
    displace!(universe, points, coords)
end


function Base.delete!(universe::Universe, point::Point)
    BasicInDelete!(universe)
    MapInDelete!(universe, point)
    NeighborInDelete!(universe, point)
    DefectAndObjectInDelete!(universe, point)
    point.debug_alive = false
end


function BasicInDelete!(universe::Universe) 
    universe.pointNum -= 1
end


function MapInDelete!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
end


function NeighborInDelete!(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(x -> x === point, neighbor.neighbors))
        PushRefreshList!(universe, neighbor.object)
    end
end


function DefectAndObjectInDelete!(universe::Universe, point::Point)
    defect = point.defect
    object = point.object
    if length(defect.pointIndexes) === 1
        deleteat!(universe.defects, findfirst(x -> x === defect, universe.defects))
        delete!(universe, object)
    else
        deleteat!(defect.pointIndexes, findfirst(x -> x === point.index, defect.pointIndexes))
        if point.type === UInt8(1)
            object.pointIndexes = defect.pointIndexes
            PushRefreshList!(universe, object)
        else
            delete!(universe, object)
        end
    end
    if point.type === UInt8(1)
        @assert(object.pointIndexes === defect.pointIndexes, "point indexes not equal")
    end
end


function displace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int32}) 
    if points[1].type === UInt8(1)
        for c in newCoords
            if c <=0 || c > universe.mapSize[1] || c > universe.mapSize[2] || c > universe.mapSize[3]
                @goto delete
            end
        end
        @goto start
        @label delete
        if rand(Uniform(0,1)) < BIAS_RATE
            for point in points
                delete!(universe, point)
            end
            return
        end
    end
    @label start
    CleanMapInDisplace!(universe, points)
    PBCCoord!(universe, newCoords)
    alivePoints = Point[]  
    BasicInDisplace!(points, newCoords)  
    for point in points
        @assert(point.debug_alive, "point not alive being displaced")
        isDeleted = OccupyInPushAndDisplace!(universe, point, point.type)
        if !isDeleted
            MapInPushAndDisplace!(universe, point)
            NeighborInDisplace!(universe, point)
            PushRefreshList!(universe, point.object)
            push!(alivePoints, point)
        else
            delete!(universe, point)
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
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(x -> x === point, neighbor.neighbors))
        PushRefreshList!(universe, neighbor.object)
    end
    point.neighbors = Point[]
    for neighborVector in NEIGHBOR_VECTORS
        coord = point.coord + neighborVector
        PBCCoord!(universe, coord)
        neighborIndex = universe.map[coord[1], coord[2], coord[3]]
        if neighborIndex > 0
            neighbor = universe.points[neighborIndex]
            push!(point.neighbors, neighbor)
            PushRefreshList!(universe, neighbor.object)
            push!(neighbor.neighbors, point)
        end
    end
end

