function RefreshObjects!(universe::Universe)
    for object in universe.objectsToRefresh
        Refresh!(object)
    end
    universe.refreshObjects = Object[]
end

function Refresh!(object::Object)
    object.isRefreshed = true
    if object.type === UInt8(2)
        RefreshVac!(object)
    else
        RefreshSIA!(object)
    end
end

function RefreshVac!(object::Object)
    vac = universe.points[object.pointIndexes[1]]
end

function VacNeighborToProbability!(neighbors::Point[])
    
