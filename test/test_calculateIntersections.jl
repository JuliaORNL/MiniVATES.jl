import MiniVATES
import MiniVATES: ScalarType, CoordType, PreallocVector
import MiniVATES: Crd4, Hist3, SquareMatrix3c
import Test: @test, @testset
import JACC

@testset "calculateIntersections" begin
    x = range(-10.0, length = 201, stop = 10.0)
    y = range(-10.0, length = 201, stop = 10.0)
    z = range(-0.1, length = 2, stop = 0.1)

    histogram = Hist3(x, y, z)

    hX = collect(x)
    kX = collect(y)
    lX = collect(z)

    doctest = MiniVATES.MDNorm(hX, kX, lX)
    maxIx = MiniVATES.maxIntersections(doctest)

    open(calc_intersections_file) do f
        line = split(strip(readline(f)), ' ')
        transform = transpose(SquareMatrix3c(parse.(ScalarType, line)))

        ndets = parse(Int, readline(f))

        intersections = PreallocVector(Vector{Crd4}(undef, maxIx))
        iPerm = PreallocVector([n for n = 1:maxIx])
        for i = 1:ndets
            line = split(strip(readline(f)), ' ')
            i_f = parse(Int, line[1])
            theta = parse(CoordType, line[2])
            phi = parse(CoordType, line[3])
            lowvalue = parse(CoordType, line[4])
            highvalue = parse(CoordType, line[5])
            if i != i_f
                i = i_f
            end
            num_intersections = parse(Int, readline(f))
            values = Vector{Vector{ScalarType}}()
            for i = 1:num_intersections
                line = split(strip(readline(f)), ' ')
                push!(values, parse.(ScalarType, line))
            end

            sortedIntersections = MiniVATES.calculateIntersections!(
                doctest,
                histogram,
                theta,
                phi,
                transform,
                lowvalue,
                highvalue,
                intersections,
                iPerm,
            )
            @test length(intersections) == num_intersections
            for i = 1:num_intersections
                for j = 1:4
                    @test isapprox(sortedIntersections[i][j], values[i][j], rtol = 0.01)
                end
            end
        end
    end
end
