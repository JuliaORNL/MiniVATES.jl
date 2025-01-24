include("test_data_constants.jl")

import MiniVATES
import MiniVATES: ScalarType, CoordType
import MiniVATES: Crd4, Hist3, SquareMatrix3c
import Test: @test, @testset
import JACC

@testset "calculateIntersections" begin
    x = range(-10.0, length = 201, stop = 10.0)
    y = range(-10.0, length = 201, stop = 10.0)
    z = range(-0.1, length = 2, stop = 0.1)

    histogram = Hist3(x, y, z)

    doctest = MiniVATES.MDNorm(x, y, z)
    maxIx = MiniVATES.maxIntersections(histogram)

    open(calc_intersections_file) do f
        line = split(strip(readline(f)), ' ')
        transform = transpose(SquareMatrix3c(parse.(ScalarType, line)))

        ndets = parse(Int, readline(f))
        doctest.intersections = MiniVATES.PreallocJaggedArray{Crd4}(ndets, maxIx)

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

            JACC.parallel_for(1,
                (i, t) -> begin
                    intersections = MiniVATES.row(t.intersections, 1)
                    MiniVATES.calculateIntersections!(
                        t.histogram,
                        t.theta,
                        t.phi,
                        t.transform,
                        t.lowvalue,
                        t.highvalue,
                        intersections,
                    )
                end,
                (
                    intersections = doctest.intersections,
                    histogram,
                    theta,
                    phi,
                    transform,
                    lowvalue,
                    highvalue
                ),
            )
            intersections = MiniVATES.row(doctest.intersections, 1)
            @test length(intersections) == num_intersections
            for i = 1:num_intersections
                for j = 1:4
                    @test isapprox(intersections[i][j], values[i][j], rtol = 0.01)
                end
            end
        end
    end
end
