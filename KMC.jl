function InputDislocationLoop(universe::Universe, pointNum::Int64, centerCoord::Vector{Float64}, directionIndex::UInt8)
    coords = HexPoints(pointNum, centerCoord, SIA_DIRECTIONS[directionIndex])
    points = Vector{Point}(undef, pointNum)
    for i in 1:size(coords)[1]
        points[i] = Point(coords[i,:])
    end
    push!(universe, points, UInt8(1), directionIndex)
end

