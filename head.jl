
const Sia_DIRECTIONS = ([Int32[1,1,1], 
                         Int32[1,1,-1], 
                         Int32[1,-1,1], 
                         Int32[1,-1,-1]])


const DISPLACE_DIRECTIONS = ([Int32[1,1,1], 
                              Int32[1,1,-1], 
                              Int32[1,-1,1], 
                              Int32[1,-1,-1],
                              Int32[-1,1,1], 
                              Int32[-1,1,-1], 
                              Int32[-1,-1,1], 
                              Int32[-1,-1,-1]])

const NEIGHBOR_VECTORS = ([Int32[2,0,0], 
                           Int32[-2,0,0], 
                           Int32[0,2,0], 
                           Int32[0,-2,0], 
                           Int32[0,0,2], 
                           Int32[0,0,-2],
                           Int32[1,1,1], 
                           Int32[1,1,-1], 
                           Int32[1,-1,1], 
                           Int32[1,-1,-1],
                           Int32[-1,1,1], 
                           Int32[-1,1,-1], 
                           Int32[-1,-1,1], 
                           Int32[-1,-1,-1]])

const NEIGHBOR_DIR = Dict([Int32[2,0,0]=>2^13, 
                     Int32[-2,0,0]=>2^12, 
                     Int32[0,2,0]=>2^11, 
                     Int32[0,-2,0]=>2^10, 
                     Int32[0,0,2]=>2^9, 
                     Int32[0,0,-2]=>2^8,
                     Int32[1,1,1]=>2^7, 
                     Int32[1,1,-1]=>2^6, 
                     Int32[1,-1,1]=>2^5, 
                     Int32[1,-1,-1]=>2^4,
                     Int32[-1,1,1]=>2^3, 
                     Int32[-1,1,-1]=>2^2, 
                     Int32[-1,-1,1]=>2^1,     
                     Int32[-1,-1,-1]=>2^0])

const DEFECT_TYPE_NAMES = ["Sia", "Vac"]

struct VacMigrationCondition
    index::UInt32
    paths::Vector{Vector{Int32}}
    startTypes::UInt8
    endTypes::Vector{UInt8}
    probabilities::Vector{Float64}
    totalProbability::Float64
end

mutable struct SiaMigrationCondition
    paths::Vector{Vector{Int32}}
    probabilities::Vector{Float64}
    totalProbability::Float64
    function SiaMigrationCondition()
        paths = Vector{Vector{Int32}}()
        probabilities = Vector{Float64}()
        totalProbability = 0
        new(paths, probabilities, totalProbability)
    end
end

const emptyVacMigrationCondition = VacMigrationCondition(
                                            UInt32(1), 
                                            Vector{Vector{Int32}}([]), 
                                            UInt8(0), UInt8[], Float64[], 0.)

mutable struct Defect
    index::UInt32    # name 
    type::UInt8
    directionIndex::UInt8
    pointIndexes::Vector{UInt32}
end


mutable struct Object
    index::UInt32
    type::UInt8
    directionIndex::UInt8
    pointIndexes::Vector{UInt32}
    isRefreshed::Bool
    vacMigrationCondition::VacMigrationCondition
    siaMigrationCondition::SiaMigrationCondition
    function Object(index::UInt32, type::UInt8, directionIndex::UInt8, pointIndexes::Vector{UInt32})
        index = index
        type = type
        directionIndex = directionIndex
        pointIndexes = pointIndexes
        emptySiaMigrationCondition = SiaMigrationCondition()
        new(index, type, directionIndex, pointIndexes, false, 
           emptyVacMigrationCondition, emptySiaMigrationCondition)
    end
end


mutable struct Point
    index::UInt32   # name and index (index nerver changed)
    coord::Vector{Int32}
    defect::Defect
    object::Object
    neighbors::Vector{Point}
    type::UInt8
    debug_alive::Bool
    function Point(coord)
        pointIndexes = UInt32[]
        defect = Defect(UInt32(0), UInt8(0), UInt8(0), pointIndexes)
        object = Object(UInt32(0), UInt8(0), UInt8(0), pointIndexes)
        neighbors = Vector{Point}[]
        new(0, coord, defect, object, neighbors, UInt8(0), true)
    end
end



mutable struct Universe
    points::Vector{Point}
    pointNum::UInt32
    defects::Vector{Defect}
    maxDefectIndex::UInt32
    objects::Vector{Object}
    maxObjectIndex::UInt32
    map::Array{UInt32,3}
    mapSize::Vector{Int32}
    nStep::Int64
    objectsToRefresh::Vector{Object}
    vacMigrationConditions::Vector{VacMigrationCondition}
    function Universe(mapSize::Vector{Int32})
        mapSize = Vector{Int32}(mapSize)
        points = Point[]
        defects = Defect[]
        objects = Object[]
        objectsToRefresh = Object[]
        map = zeros(UInt32, mapSize[1], mapSize[2], mapSize[3])
        vacMigrationConditions = VacMigrationCondition[]
        new(points, UInt32(0), defects, UInt32(0), objects, UInt32(0), map, mapSize, 0, objectsToRefresh, vacMigrationConditions)
    end
end



function Base.display(point::Point)
    print("$(point.index) $(DEFECT_TYPE_NAMES[point.type]) $(point.defect.index) $(point.object.index) \
          ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
    println()
end


function Base.display(points::Vector{Point})
    println("$(length(points))-point array:")
    println("id type defect object (x y z)")  # id is index
    for point in points
        display(point)
    end
    println()
end



function Base.display(universe::Universe, defect::Defect)
    direction = Sia_DIRECTIONS[defect.directionIndex]
    print("id: $(defect.index)  type: $(DEFECT_TYPE_NAMES[defect.type]) ")
    if defect.type === UInt8(1)
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

function Base.display(c::VacMigrationCondition)
    print("id: $(c.index) | probabilities: $(c.probabilities) | total probability: $(c.totalProbability)")
end