using StatsBase
using Crayons
include("geometry.jl")
using .Geometry

const directionDir = Int64[[1,1,1],
                      [1,1,-1],
                      [1,-1,1],
                      [1,-1,-1],
                      [0,0,0]]

                      
mutable struct Point 
    id::Int64
    index::UInt32
    type::UInt8
    coord::Vector{Int64}
    neighbors::Vector{Point}
    directionIndex::Vector{Int64}
    defectIndex::Int64
    function Point(type::UInt8, coord::Vector{Int64}, directionIndex::Vector{Int64})
        @assert(type==1 || type==2, "Bad point type")
        @assert((type==2 && directionIndex==5) || type==1 && directionIndex<5, "Bad direction index")
        neighbors = Point[]
        new(0, UInt32(0), type, coord, neighbors, directionIndex, 0)
    end
end


const pointTypeNames = ["SIA", "Vac"]
function Base.display(point::Point)
    print("$(point.id) $(point.index) $(pointTypeNames[point.type]) \
          ($(point.coord[1]) $(point.coord[2]) $(point.coord[3]))")
end


function Base.display(points::Vector{Point})
    println("$(length(points))-point universe:")
    println("id index type (x y z)")
    for point in points
        display(point)
        println()
    end
end


mutable struct Defect
    index::UInt32
    id :: Int64
    type::UInt8
    points::Vector{Point}
    directionIndex::Vector{Int64}
    function Defect()
        points = Point[]
        new(UInt32(0), 0, UInt8(0), points, [0,0,0])
    end
end


mutable struct Universe
    mapsize::Vector{Int64}
    map::Array{UInt32, 3}
    points::Vector{Point}
    pointNum::UInt32
    maxPointID::Int64
    defects::Vector{Defect}
    defectNum::UInt32
    maxDefectID::Int64
    nstep::Int64
    function Universe(mapsize::Vector{Int64})
        map = zeros(UInt32, mapsize[1], mapsize[2], mapsize[3]) # N x half a lattice constant
        points = Point[]
        defects = Defect[]
        new(mapsize, map, points, UInt32(0), 0, defects, UInt32(0), 0, 0)
    end
end


# dump the map in lammps dump format
function Dump(universe::Universe, filename::String, mode::String)
    file = open(filename, mode)
    write(file, "ITEM: TIMESTEP\n")
    write(file, "$(universe.nstep)\n")
    write(file, "ITEM: NUMBER OF ATOMS\n")
    write(file, "$(universe.pointNum)\n")
    write(file, "ITEM: BOX BOUNDS pp pp pp\n")
    write(file, "1 $(universe.mapsize[1])\n")
    write(file, "1 $(universe.mapsize[2])\n")
    write(file, "1 $(universe.mapsize[3])\n")
    write(file, "ITEM: ATOMS id type x y z\n")
    for i=1:universe.pointNum
        point = universe.points[i]
        write(file, "$(point.id) $(point.type) $(point.coord[1]) $(point.coord[2]) $(point.coord[3])\n")
    end
    close(file)
end


function Base.display(universe::Universe, zorder::Int64)
    map2D = universe.map[:,:,zorder]
    print("  ")
    for i=1:universe.mapsize[1]
        print(" $(i%10)")
    end
    println()
    for i=1:universe.mapsize[1]
        print(" $(i%10)")
        for j=1:universe.mapsize[2]
            if isodd(zorder)
                if isodd(i) && isodd(j)
                    c = map2D[i,j] > 0 ? "$(map2D[i,j])" : " "
                    print(Crayon(background=:dark_gray, foreground=:white, bold=true), " ", c)
                    print(Crayon(reset=true),"")
                else
                    print(Crayon(background=:default, bold=true),"  ")
                    print(Crayon(reset=true),"")
                end
            else 
                if iseven(i) && iseven(j)
                    c = map2D[i,j] > 0 ? "$(map2D[i,j])" : " "
                    print(Crayon(background=:dark_gray, foreground=:white, bold=true), " ", c)
                    print(Crayon(reset=true),"")
                else
                    print(Crayon(background=:default, bold=true),"  ")
                    print(Crayon(reset=true),"")
                end
            end
        end
        println()
    end
end


function RefreshFile(filename::String)
    file = open(filename, "w")
    close(file)
end


function Base.push!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] == 0, "points overlap!")
    PushBasic!(universe, point)
    PushMap!(universe, point)
    PushNeighbors!(universe, point)
    PushSVReact!(universe, point)
end

function PushBasic!(universe::Universe, point::Point)
    universe.pointNum += 1
    universe.maxPointID += 1
    point.index = universe.pointNum
    point.id = universe.maxPointID
    push!(universe.points, point)
end

function PushMap!(universe::Universe, point::Point)
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = point.index
end

function PushNeighbors!(universe::Universe, point::Point)
    for i in [-1,1]
        for j in [-1,1]
            for k in [-1,1]
                coord = point.coord + [i,j,k]
                try
                    if universe.map[coord[1], coord[2], coord[3]] > 0
                        neighbor = universe.points[universe.map[coord[1], coord[2], coord[3]]]
                        push!(point.neighbors, neighbor)
                        push!(neighbor.neighbors, point)
                    end
                catch BoundsError
                    continue
                end
            end
        end
    end
end


function PushSVReact!(universe::Universe, point::Point)
    reversePoints = Point[]
    for neighbor in point.neighbors
        if neighbor.type != point.type
            push!(reversePoints, neighbor)
        end
    end
    if !isempty(reversePoints)
        neighbor = sample(reversePoints)
        delete!(universe, neighbor)
        delete!(universe, point)
    end
end


function PushDefect!(universe::Universe, point::Point)  #Defecting after a push
    if point.type == 0 || (point.type == 1 && isempty(point.neighbors))
        defect = Defect()
        defect.points = [point]
        defect.directionIndex = point.directionIndex
        defect.type = point.type
        push!(universe, defect)
    else
        defectIndex = point.neighbors[1].defectIndex
        defect = universe.defects[defectIndex]
        push!(defect.points, point)
        mergedIndexes = [defectIndex]
        for neighbor in point.neighbors[2:end]
            if !(neighbor.defectIndex in mergedIndexes) 
                push!(mergedIndexes, defectIndex)
                Merge!(universe::Universe, defect, universe.defects[neighbor.defectIndex])
                push!(neighbor.defectIndex, mergedIndexes)
            end
        end
        Rearrange!(universe, defect)
    end
end

function Rearrange!(universe::Universe, defect::Defect)
    directionIndexList = [0,0,0,0,0]
    for p in defect.points
        directionIndexList[p.directionIndex] += 1
    end
    maxDirectionIndex = argmax(directionIndexList)
    defect.directionIndex = maxDirectionIndex
    aveCoord = [0,0,0]
    for p in defect.points
        p.directionIndex = maxDirectionIndex
        aveCoord[1] += p.coord[1]
        aveCoord[2] += p.coord[2]
        aveCoord[3] += p.coord[3]
    end
    aveCoord[1] /= UInt32(round(length(defect.points)))
    aveCoord[2] /= UInt32(round(length(defect.points)))
    aveCoord[3] /= UInt32(round(length(defect.points)))
    coords = HexPoints(length(defect.points), aveCoord, directionDir[maxDirectionIndex])
    Displace!(universe, defect.points, coords)
end


function Base.push!(defect::Defect, point::Point)
    push!(defect.points, point)
    point.defectIndex = defect.index
end
        
function Merge!(universe::Universe, defect1::Defect, defect2::Defect)
    defect1.points = vcat(defect1.points, defect2.points)
    for point in defect2.points
        point.defectIndex = defect1.index
    end
    delete!(universe, defect2)
end

function Base.delete!(universe::Universe, defect::Defect)
    deleteat!(universe.defects, defect.index)
    for defect in universe.defects[defect.index:end]
        SetDefectIndex!(defect, defect.index - 1)
    end
end

function SetDefectIndex(universe::Universe, defect::Defect, index::Int64)
    defect.index = index
    for p in defect.points
        p.defectIndex = index
    end
end


function Base.push!(universe::Universe, defect::Defect)
    push!(universe.defects, defect)
    universe.defectNum += 1
    defect.index = universe.defectNum
    universe.maxDefectID += 1
    defect.id = universe.maxDefectID
    for point in defect.points
        point.defectIndex = defect.index
    end
end



function Base.delete!(universe::Universe, point::Point)
    DeleteNeighbors!(universe, point)
    DeleteBasic!(universe, point)
    DeleteMap!(universe, point)
    DeleteDefect!(universe, point)
end

function DeleteBasic!(universe::Universe, point::Point)
    universe.pointNum -= 1
    deleteat!(universe.points, point.index)
    for point in universe.points[point.index:end]
        point.index -= 1
    end
end

function DeleteDefect!(universe::Universe, point::Point)
    defect = universe.defects[point.defectIndex]
    if length(defect.points) == 1
        delete!(universe, defect)
    else
        deleteat!(defect.points, indexof(defect.points, point))
    end
end

function DeleteMap!(universe::Universe, point::Point)
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] > 0, "point not in map!")
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = 0
    for p in universe.points[point.index:end]
        universe.map[p.coord[1], p.coord[2], p.coord[3]] -= 1
    end
end


function DeleteNeighbors!(universe::Universe, point::Point)
    for neighbor in point.neighbors
        deleteat!(neighbor.neighbors, findfirst(p->p.id==point.id,neighbor.neighbors))
    end
end


function Displace!(universe::Universe, point::Point, newCoord::Vector{Int64})
    @assert(universe.map[point.coord[1], point.coord[2], point.coord[3]] > 0, "point not in map!")
    @assert(universe.map[newCoord[1], newCoord[2], newCoord[3]] == 0, "points overlap!")
    DisplaceBasic!(universe, point, newCoord)
    DisplaceMap!(universe, point, newCoord)
end

function DisplaceBasic!(universe::Universe, point::Point, newCoord::Vector{Int64})
    point.coord = newCoord
end

function DisplaceMap!(universe::Universe, point::Point, newCoord::Vector{Int64})
    universe.map[point.coord[1], point.coord[2], point.coord[3]] = UInt32(0)
    point.coord = newCoord
    universe.map[newCoord[1], newCoord[2], newCoord[3]] = point.index
end

function DisplaceNeighbors!(universe::Universe, point::Point) 
    # use it after a lot of points have been displaced
    neighborIndexes = UInt32[]
    neighbors = Point[]
    neighborsBefore = point.neighbors
    for i in [-1,1]
        for j in [-1,1]
            for k in [-1,1]
                coord = point.coord + [i,j,k]
                try
                    index = universe.map[coord[1], coord[2], coord[3]] 
                catch BoundsError
                    continue
                end
                if index > 0
                    neighbor = universe.points[index]
                    push!(neighborIndexes, index)
                    push!(neighbors, neighbor)
                end
            end
        end
    end
    neighborIndexesBefore = [p.index for p in neighborsBefore]
    dropedNeighborIndexes = setdiff(neighborIndexesBefore, neighborIndexes)
    newNeighborIndexes = setdiff(neighborIndexes, neighborIndexesBefore)
    for index in dropedNeighborIndexes
        neighbor = universe.points[index]
        deleteat!(neighbor.neighbors, findfirst(p->p.id==point.id,neighbor.neighbors))
    end
    for index in newNeighborIndexes
        neighbor = universe.points[index]
        push!(neighbor.neighbors, point)
    end
    point.neighbors = neighbors
end

function Displace!(universe::Universe, points::Vector{Point}, newCoords::Matrix{Int64})
    for i in 1:length(points)
        Displace!(universe, point, newCoords[i,:])
    end
    DisplaceNeighbors!(universe, points)
end

function DisplaceSVReact!(universe::Universe, points::Vector{Point})
    for point in points
        for neighbor in point.neighbors
            reversePoints = Point[]
            if neighbor.type != point.type
                push!(reversePoints, neighbor)
            end
            if !isempty(reversePoints)
                neighbor = sample(reversePoints)
                delete!(universe, neighbor)
                delete!(universe, point)
            end
        end
    end
end

function DisplaceDefect!(universe::Universe, points::Vector{Point})
            


universe = Universe([10,10,10])
point1 = Point(UInt8(1), [3,3,3], [1,1,1])
point2 = Point(UInt8(1), [4,4,4], [1,1,1])
point3 = Point(UInt8(1), [2,2,2], [1,1,1])


filename = "/mnt/c/Users/xuke/Desktop/test.dump"
RefreshFile(filename)
push!(universe, point1)
push!(universe, point2)
push!(universe, point3)
Dump(universe, filename, "a")

