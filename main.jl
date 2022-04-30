using StatsBase
using Random
using Distributions
Random.seed!(1234)

include("geometry.jl")
using .Geometry
include("reaction.jl")
include("KMC.jl")

function test1(universe::Universe)
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


function test2(universe::Universe)
    for i in 1:1000
        universe.nStep += 1
        println("step: ", universe.nStep)
        point = Point(Vector{Int32}([150,150,150]))
        push!(universe, [point], UInt8(2), UInt8(1))
        if universe.nStep % 10 == 0
            Dump(universe, fileName, "a")
        end
    end
    universe
end


function test3(universe::Universe)
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
end

function test4(universe::Universe)
    InputDislocationLoop(universe, 100, Vector{Float64}([150,150,150]), UInt8(2))
    Dump(universe, filename, "a")
end


mapSize = Vector{UInt32}([300,300,300])
universe = Universe(mapSize)
fileName = "/mnt/c/Users/xuke/Desktop/test3.dump"
RefreshFile(fileName)

test3(universe)