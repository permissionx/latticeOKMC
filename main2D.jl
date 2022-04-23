using StatsBase
mutable struct Point
    index::UInt32
    type::Int8 # 1 for sia, 2 for vac
    coord::Vector{Int64}
    neighbors::Vector{Point}
    groupIndex::UInt32
    function Point(type::Int64, coord::Vector{Int64})
        type = Int8(type)
        @assert(type==1 || type==2, "Bad point type")
        neighbors = Point[]
        new(0, type, coord, neighbors)
    end
end

const pointTypeNames = ["SIA", "Vac"]
function Base.show(io::IO, point::Point)
    print(io, "$(pointTypeNames[point.type]) ($(point.coord[1]) $(point.coord[2]))")
end

mutable struct Group 
    index::UInt32
    type::UInt8
    pointContainer::Vector{Point}
    function Group()
        new(0,0,Point[])
    end
end


function SetGroupIndex!(group::Group, index::UInt32)
    group.index = index
    for point in group.pointContainer
        point.groupIndex = index
    end
end

    

const MAPSIZE = 5
mutable struct Universe
    map::Array{UInt32,2}
    pointContainer::Vector{Point}
    pointNum::UInt32
    groups::Vector{Group}
    groupNum::UInt32
    function Universe()
        map = zeros(UInt32,MAPSIZE,MAPSIZE)
        pointContainer = Point[]
        groups = Group[]
        new(map, pointContainer, UInt32(0), groups, UInt32(0))
    end
end

function Group_Add!(universe::Universe, point::Point)
    neighbors = point.neighbors
    if isempty(neighbors)
        universe.groupNum += 1
        group = Group()
        group.pointContainer = [point]
        group.type = point.type
        SetGroupIndex!(group, universe.groupNum)
        push!(universe.groups, group)
    else
        gIndex1 = neighbors[1].groupIndex
        push!(universe.groups[gIndex1].pointContainer, point)
        point.groupIndex = gIndex1
        @assert(universe.groups[gIndex1].type == point.type, "point type and group type not consistent")
        for neighbor in neighbors[2:end]
            gIndexLatter = neighbor.groupIndex
            if gIndexLatter != gIndex1
                @assert(universe.groups[gIndex1].type == universe.groups[gIndexLatter].type,
                "merge groups with different types")
                universe.groups[gIndex1] = vcat(universe.groups[gIndex1], universe.groups[gIndexLatter])
                SetGroupIndex!(universe.groups[gIndex1], gIndex1)
                DeleteGroup!(universe, gIndexLatter)
            end
        end
    end
end

function Group_Delete!(universe::Universe, point::Point)
    neighbors = point.neighbors
    labels = zeros()



function DeleteGroup!(universe::Universe, index::UInt32)
    deleteat!(universe.groups, index)
    for group in universe.groups[index+1, end]
        SetGroupIndex!(group, group.index-1)
    end
end


function SetMapVaule!(universe::Universe, coord::Vector{Int64}, index::UInt32)
    @assert(!(index!=0 && universe.map[coord[1], coord[2]]!=0), "Points overlaped! $index $(universe.map[coord[1], coord[2]])")
    universe.map[coord[1], coord[2]] = index
end


function Base.push!(universe::Universe, point::Point)
    universe.pointNum += 1
    push!(universe.pointContainer, point)
    SetMapVaule!(universe, point.coord, universe.pointNum)
    point.index = universe.pointNum
    neighbors = FindNeighbor(universe, point)
    SetNeighbors!(point, neighbors)
    for neighbor in neighbors
        PushNeighbor!(neighbor, point)
    end
    React_AddAndRelocate!(universe, point)
    Group_Add!(universe, point)
end

function Base.delete!(universe::Universe, point::Point)
    neighbors = FindNeighbor(universe, point)
    for neighbor in neighbors
        DeleteNeighbor!(neighbor, point)
    end
    SetMapVaule!(universe, point.coord, UInt32(0))
    index = point.index
    for p in universe.pointContainer[index+1:end]
        universe.map[p.coord[1], p.coord[2]] -= 1
        p.index -= 1
    end
    deleteat!(universe.pointContainer, index)
    universe.pointNum -= 1
end

function Relocate!(universe::Universe, point::Point, newCoord::Vector{Int64})
    SetMapVaule!(universe, point.coord, UInt32(0))
    point.coord = newCoord
    SetMapVaule!(universe, newCoord, point.index)
    neighbors = point.neighbors
    for neighbor in neighbors
        DeleteNeighbor!(neighbor, point)
    end
    neighbors = FindNeighbor(universe, point)
    SetNeighbors!(point, neighbors)
    for neighbor in neighbors
        PushNeighbor!(neighbor, point)
    end
    React_AddAndRelocate!(universe, point)
end

function FindNeighbor(universe::Universe, point::Point)
    coord = point.coord
    neighbors = Point[]
    for x in [-1, 0, 1]
        for y in [-1, 0, 1]
            if x == 0 && y == 0 
                continue
            end
                neighborCoord = [coord[1] + x, coord[2] + y]
            try
                index = universe.map[neighborCoord[1], neighborCoord[2]]
                if index > 0 
                    push!(neighbors, universe.pointContainer[index])
                end
            catch BoundsError
                continue
            end
        end
    end
    neighbors
end

function SetNeighbors!(point::Point, neighbors::Vector{Point})
    point.neighbors = neighbors
end

function PushNeighbor!(point::Point, neighbor::Point)
    push!(point.neighbors, neighbor)
end

function DeleteNeighbor!(point::Point, neighbor::Point)
    deleteat!(point.neighbors, findall(x->x==neighbor, point.neighbors))
end

function React_AddAndRelocate!(universe::Universe, point::Point)
    # to do: delete reverse point by random
    reversePoints = Point[]
    for neighbor in point.neighbors
        if neighbor.type + point.type == 3
            push!(reversePoints, neighbor)
        end
    end
    if !isempty(reversePoints)
        neighbor = sample(reversePoints)
        delete!(universe, neighbor)
        delete!(universe, point)
    end
end


function SpeedTest()
    universe = Universe()
    point1 = Point(1, [2,3])
    push!(universe, point1)
    point2 = Point(1, [2,4])
    push!(universe, point2)
    point3 = Point(1, [2,5])
    push!(universe, point3)
end

#SpeedTest()
"
vac = Point(2,[4,4])
push!(universe, vac)
Relocate!(universe, vac, [3,4])
map = universe.map
"



