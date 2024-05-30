import MiniVATES
import MiniVATES: SizeType, ScalarType, CoordType
import MiniVATES: SquareMatrix3c, C3, Crd3, Crd4
import MiniVATES: Hist3, atomic_push!, binweights, reset!, maxIntersections
import MiniVATES: PreallocVector, Array1
import Test: @test, @testset
import HDF5
import JACC

import Profile
import ProfileCanvas
import PProf
using StatProfilerHTML

# @testset "calculateIntersections" begin
#     x = range(-10.0, length = 201, stop = 10.0)
#     y = range(-10.0, length = 201, stop = 10.0)
#     z = range(-0.1, length = 2, stop = 0.1)

#     histogram = Hist3(x, y, z)

#     hX = collect(x)
#     kX = collect(y)
#     lX = collect(z)

#     doctest = MiniVATES.MDNorm(hX, kX, lX)
#     maxIx = maxIntersections(doctest)

#     open(calc_intersections_file) do f
#         line = split(strip(readline(f)), ' ')
#         transform = transpose(SquareMatrix3c(parse.(ScalarType, line)))

#         ndets = parse(Int, readline(f))

#         intersections = PreallocVector(Vector{Crd4}(undef, maxIx))
#         iPerm = PreallocVector([n for n = 1:maxIx])
#         for i = 1:ndets
#             line = split(strip(readline(f)), ' ')
#             i_f = parse(Int, line[1])
#             theta = parse(CoordType, line[2])
#             phi = parse(CoordType, line[3])
#             lowvalue = parse(CoordType, line[4])
#             highvalue = parse(CoordType, line[5])
#             if i != i_f
#                 i = i_f
#             end
#             num_intersections = parse(Int, readline(f))
#             values = Vector{Vector{ScalarType}}()
#             for i = 1:num_intersections
#                 line = split(strip(readline(f)), ' ')
#                 push!(values, parse.(ScalarType, line))
#             end

#             sortedIntersections = MiniVATES.calculateIntersections!(
#                 doctest,
#                 histogram,
#                 theta,
#                 phi,
#                 transform,
#                 lowvalue,
#                 highvalue,
#                 intersections,
#                 iPerm,
#             )
#             @test length(intersections) == num_intersections
#             for i = 1:num_intersections
#                 for j = 1:4
#                     @test isapprox(sortedIntersections[i][j], values[i][j], rtol = 0.01)
#                 end
#             end
#         end
#     end
# end

function kernel1_1D(i, t)
    @inbounds begin
        if !t.use_dets[i]
            return
        end

        detID = t.detIDs[i]
        wsIdx = get(t.fluxDetToIdx, detID, nothing)
        if wsIdx == nothing
            return
        end

        sortedIntersections = MiniVATES.calculateIntersections!(
            t.doctest,
            t.signal,
            t.thetaValues[i],
            t.phiValues[i],
            t.transform,
            t.lowValues[i],
            t.highValues[i],
            t.intersections[i],
            t.iPerm[i],
        )

        if isempty(sortedIntersections)
            return
        end

        MiniVATES.calculateDiffractionIntersectionIntegral!(
            sortedIntersections,
            t.integrFlux_x,
            t.integrFlux_y[wsIdx],
            t.yValues[i],
        )

        solid::ScalarType = t.protonCharge
        if t.haveSA
            saIdx = t.solidAngDetToIdx[detID]
            saFactor = t.solidAngleWS[saIdx][1]
            solid *= saFactor
        end
        MiniVATES.calculateSingleDetectorNorm!(
            t.doctest,
            sortedIntersections,
            solid,
            t.yValues[i],
            t.signal,
        )
    end
end

@testset "calculateDiffractionIntersectionIntegral" begin
    x = range(-10.0, length = 201, stop = 10.0)
    y = range(-10.0, length = 201, stop = 10.0)
    z = range(-0.1, length = 2, stop = 0.1)

    signal = Hist3(x, y, z)

    hX = collect(x)
    kX = collect(y)
    lX = collect(z)

    doctest = MiniVATES.MDNorm(hX, kX, lX)

    # rot file
    rotData = MiniVATES.loadRotationData(rot_nxs_file)
    transforms = MiniVATES.makeTransforms(rotData)

    # sa file
    saData = MiniVATES.loadSAData(sa_nxs_file)

    haveSA = true

    # flux file
    fluxData = MiniVATES.loadFluxData(flux_nxs_file)
    ndets = fluxData.ndets

    # use dets file
    use_dets = parse.(Bool, readlines(use_dets_file))
    for v in use_dets
        @test v == true
    end

    # event file
    eventsData = MiniVATES.loadEventsData(event_nxs_file)

    maxIx = maxIntersections(doctest)
    intersections = [PreallocVector(Vector{Crd4}(undef, maxIx)) for i = 1:ndets]
    yValues = [PreallocVector(Vector{ScalarType}(undef, maxIx)) for i = 1:ndets]
    iPerm = [PreallocVector([n for n = 1:maxIx]) for i = 1:ndets]

    function launch_kernel1()
        @time begin
            for n = 1:length(transforms)
                JACC.parallel_for(
                    ndets,
                    kernel1_1D,
                    (
                        transform = transforms[n],
                        use_dets = use_dets,
                        intersections,
                        iPerm,
                        yValues,
                        eventsData.detIDs,
                        fluxData.fluxDetToIdx,
                        doctest,
                        signal,
                        eventsData.thetaValues,
                        eventsData.phiValues,
                        eventsData.lowValues,
                        eventsData.highValues,
                        fluxData.integrFlux_x,
                        fluxData.integrFlux_y,
                        saData.solidAngDetToIdx,
                        saData.solidAngleWS,
                        eventsData.protonCharge,
                        haveSA,
                    ),
                )
            end
        end
    end
    launch_kernel1()
    reset!(signal)
    launch_kernel1()
    # Profile.clear()
    # Profile.@profile launch_kernel1()
    # statprofilehtml()
    # return nothing

    # norm file
    file = HDF5.h5open(norm_nxs_file, "r")
    normGroup = file["MDHistoWorkspace"]["data"]
    data = read(normGroup["signal"])
    dims = size(data)
    @test length(dims) == 3
    @test dims[1] == 200
    @test dims[2] == 200
    @test dims[3] == 1
    close(file)

    data2d = data[:, :, 1]
    max_signal = maximum(signal.weights)
    # max_signal = 0.0
    ref_max = 0.0
    for i = 1:dims[2]
        for j = 1:dims[1]
            @test isapprox(data2d[i, j], binweights(signal)[i, j, 1], atol = 1.0e+6)
            ref_max = max(ref_max, data2d[j, i])
            # max_signal = max(max_signal, binweights(signal)[i, j, 1])
        end
    end
    @test isapprox(max_signal, ref_max, atol = 1.0e+5)

    h = Hist3(x, y, z)

    # TODO:
    # - turn this into a package function
    # - write loader functions for the various hdf5 files
    # - make test to run this kernel on CUDA
    # - introduce MPI

    @time MiniVATES.binEvents!(h, eventsData.events, transforms)
    reset!(h)
    @time MiniVATES.binEvents!(h, eventsData.events, transforms)
    # Profile.clear()
    # Profile.@profile launch_kernel2()
    # statprofilehtml()

    fio = open("meow.txt", "w")
    for j = 1:dims[2]
        for i = 1:dims[1]
            println(fio, binweights(h)[i, j, 1] / binweights(signal)[i, j, 1])
        end
    end
end
