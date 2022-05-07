function RefreshObjects!(universe::Universe)
    for object in universe.objectsToRefresh
        Refresh!(object)
    end
    universe.objectsToRefresh = Object[]
end

function Refresh!(object::Object)
    object.isRefreshed = true
    if object.type === UInt8(2)
        RefreshVac!(object)
    else
        RefreshSia!(object)
    end
end

function VacNeighborConditionIndex(object::Object)
    point = universe.points[object.pointIndexes[1]]
    index = UInt32(0)
    for neighbor in point.neighbors 
        try
            index += NEIGHBOR_DIR[neighbor.coord - point.coord] 
        catch KeyError
            delta = neighbor.coord - point.coord
            for i in 1:3
                if delta[i] > Int32(2)
                    delta[i] -= mapSize[i]
                elseif delta[i] < Int32(-2)
                    delta[i] += mapSize[i]
                end
            end
            index += NEIGHBOR_DIR[delta] 
        end
    end
    index+UInt32(1)
end

function RefreshVac!(object::Object)
    migrationConditionIndex = VacNeighborConditionIndex(object)
    object.vacMigrationCondition = universe.vacMigrationConditions[migrationConditionIndex]
end


function RefreshSia!(object::Object) 
    nSia = length(object.pointIndexes)
    path = Sia_DIRECTIONS[object.directionIndex]
    paths = [path, -path]
    p = 1
    probabilities = [p, p]
    totalProbability = p*2
    object.siaMigrationCondition.paths = paths
    object.siaMigrationCondition.probabilities = probabilities
    object.siaMigrationCondition.totalProbability = totalProbability
end


function IterStep!(universe::Universe)
    RefreshObjects!(universe)
    object = RandomAnObject(universe)
    path = RandomPath(object)
    displace!(universe, object, path)
end

function RandomAnObject(universe::Universe)  # todo: should be optimized
    probabilityList = Float64[]
    for object in universe.objects
        if object.type === UInt8(2)
            push!(probabilityList, object.vacMigrationCondition.totalProbability)
        else
            push!(probabilityList, object.siaMigrationCondition.totalProbability)
        end
    end
    weights = Weights(probabilityList)
    object = sample(universe.objects, weights)
    object
end

function RandomPath(object::Object)
    if object.type === UInt8(2)
        weights = Weights(object.vacMigrationCondition.probabilities)
        path = sample(object.vacMigrationCondition.paths, weights)
    else
        weights = Weights(object.siaMigrationCondition.probabilities)
        path = sample(object.siaMigrationCondition.paths, weights)
    end
    path
end

function displace!(universe::Universe, object::Object, path::Vector{Int32})
    nPoints = length(object.pointIndexes)
    newCoords = Matrix{Int32}(undef, nPoints, 3)
    points = Point[]
    for i in 1:nPoints
        point = universe.points[object.pointIndexes[i]]
        push!(points, point)
        newCoords[i,:] = point.coord + path
    end
    displace!(universe, points, newCoords) 
end

