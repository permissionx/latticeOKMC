module VacMigration

function GenerateNeighborConditions()
    neighborConditions = Vector{Vector{Vector{Int32}}}([])
    neighbors = Vector{Vector{Int32}}([])
    IterNeighbor!(neighborConditions, neighbors,1)
    neighborConditions
end

function IterNeighbor!(neighborConditions::Vector{Vector{Vector{Int32}}}, neighbors::Vector{Vector{Int32}}, n::Int64)
    if n > 14
        push!(neighborConditions, neighbors)
        return
    else
        IterNeighbor!(neighborConditions, neighbors, n+1)
        neighbor = Main.NEIGHBOR_VECTORS[n]
        IterNeighbor!(neighborConditions, vcat(neighbors, [neighbor]), n+1)
    end
end

function Paths(neighbors::Vector{Vector{Int32}})
    paths = setdiff(Main.NEIGHBOR_VECTORS, neighbors)
    types = UInt8[]
    for path in paths
        for neighbor in neighbors
            delta = abs.(path-neighbor)
            if delta == Int32[1,1,1]
                push!(types, UInt8(1))
                @goto finish
            end
        end
        for neighbor in neighbors
            delta = abs.(path-neighbor)
            if delta == Int32[2,0,0] || delta == Int32[0,2,0] || delta == Int32[0,0,2]
                push!(types, UInt8(2))
                @goto finish
            end
        end
        push!(types, UInt8(3))
        @label finish
    end
    paths, types
end

function MigrationConditionIndex(neighbors::Vector{Vector{Int32}})
    index = Int32(0)
    for neighbor in neighbors
        index += Main.NEIGHBOR_DIR[neighbor] 
    end
    index+UInt32(1)
end    

function StartType(neighbors::Vector{Vector{Int32}})
    for neighbor in neighbors
        absNeighbor = abs.(neighbor)
        if absNeighbor == Int32[1,1,1]
            return UInt8(1)
        end
    end
    for neighbor in neighbors
        absNeighbor = abs.(neighbor)
        if absNeighbor == Int32[2,0,0] || absNeighbor == Int32[0,2,0] || absNeighbor == Int32[0,0,2]
            return UInt8(2)
        end
    end
    UInt8(3)
end 

function Probabilities(endTypes::Vector{UInt8}, startType::UInt8)
    ps = Float64[]
    for endType in endTypes
        if startType == UInt8(3)
            barrier = 1
        elseif startType == UInt8(1)
            if endType == UInt8(1)
                barrier = 100000
            elseif endType == UInt8(2)
                barrier = 100000
            else
                barrier = 100000
            end
        else  # startType == UInt8(2)
            if endType == UInt8(3)
                barrier = 100000
            elseif endType == UInt8(2)
                barrier = 100000
            else
                barrier = 1
            end
        end
        p = 0.1/barrier
        push!(ps, p)
    end
    totalProbability = sum(ps)
    ps, totalProbability
end


function MigrationCondition(neighbors::Vector{Vector{Int32}})
    paths, endTypes = Paths(neighbors)
    index = MigrationConditionIndex(neighbors)
    startType = StartType(neighbors)
    probabilities, totalProbability = Probabilities(endTypes, startType)
    migrationCondition = Main.VacMigrationCondition(index, paths, startType, endTypes, probabilities, totalProbability)
    migrationCondition
end

function GenerateConditions!(universe::Main.Universe)
    neighborConditions = GenerateNeighborConditions()
    conditions = Main.VacMigrationCondition[]
    for neighbors in neighborConditions
        condition = MigrationCondition(neighbors)
        push!(conditions, condition)
    end
    universe.vacMigrationConditions = conditions
end
    
end

module Init
export Initialization!
function Initialization!(universe::Main.Universe)
    Main.VacMigration.GenerateConditions!(universe)
    println("Initialization completed!")
end
end

