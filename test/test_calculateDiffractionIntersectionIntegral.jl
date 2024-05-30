import MiniVATES
import MiniVATES: ScalarType, CoordType
import MiniVATES: SquareMatrix3c, Crd4
import MiniVATES: Hist3, binweights, reset!
import MiniVATES: PreallocVector
import Test: @test, @testset
import HDF5
import JACC

import Profile
import StatProfilerHTML: statprofilehtml

function kernel1_1D(i, t)
    @inbounds begin
        if t.skip_dets[i]
            return nothing
        end

        detID = t.detIDs[i]
        wsIdx = get(t.fluxDetToIdx, detID, nothing)
        if wsIdx == nothing
            return nothing
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
            return nothing
        end

        MiniVATES.calculateDiffractionIntersectionIntegral!(
            sortedIntersections,
            t.integrFlux_x,
            t.integrFlux_y[wsIdx],
            t.yValues[i],
        )

        saIdx = t.solidAngDetToIdx[detID]
        saFactor = t.solidAngleWS[saIdx][1]
        solid::ScalarType = t.protonCharge * saFactor

        MiniVATES.calculateSingleDetectorNorm!(
            t.doctest,
            sortedIntersections,
            solid,
            t.yValues[i],
            t.signal,
        )

        return nothing
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
    exData = MiniVATES.loadExtrasData(rot_nxs_file)
    m_W = SquareMatrix3c([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])
    transforms = MiniVATES.makeRotationTransforms(exData, m_W)
    transforms2 = MiniVATES.makeTransforms(exData, m_W)

    # sa file
    saData = MiniVATES.loadSolidAngleData(sa_nxs_file)

    # flux file
    fluxData = MiniVATES.loadFluxData(flux_nxs_file)
    ndets = fluxData.ndets

    # event file
    eventData = MiniVATES.loadEventData(event_nxs_file)

    maxIx = MiniVATES.maxIntersections(doctest)
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
                        skip_dets = exData.skip_dets,
                        intersections,
                        iPerm,
                        yValues,
                        eventData.detIDs,
                        fluxData.fluxDetToIdx,
                        doctest,
                        signal,
                        eventData.thetaValues,
                        eventData.phiValues,
                        eventData.lowValues,
                        eventData.highValues,
                        fluxData.integrFlux_x,
                        fluxData.integrFlux_y,
                        saData.solidAngDetToIdx,
                        saData.solidAngleWS,
                        eventData.protonCharge,
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
    ref_max = 0.0
    for i = 1:dims[2]
        for j = 1:dims[1]
            @test isapprox(data2d[i, j], binweights(signal)[i, j, 1], atol = 1.0e+6)
            ref_max = max(ref_max, data2d[j, i])
        end
    end
    @test isapprox(max_signal, ref_max, atol = 1.0e+5)

    h = Hist3(x, y, z)

    # TODO:
    # - introduce MPI
    # - make test to run this kernel on CUDA

    @time MiniVATES.binEvents!(h, eventData.events, transforms2)
    reset!(h)
    @time MiniVATES.binEvents!(h, eventData.events, transforms2)
    # Profile.clear()
    # Profile.@profile launch_kernel2()
    # statprofilehtml()

    open("meow.txt", "w") do fio
        for j = 1:dims[2]
            for i = 1:dims[1]
                println(fio, binweights(h)[i, j, 1] / binweights(signal)[i, j, 1])
            end
        end
    end
end
