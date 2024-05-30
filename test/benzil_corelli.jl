include("test_data_constants.jl")

import MiniVATES
import MiniVATES: ScalarType, CoordType
import MiniVATES: Hist3, SquareMatrix3c, Crd4, PreallocVector
import Test: @testset
import JACC

#TODO: make this a package function
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

@testset "benzil_corelli" begin
    x = range(start = -7.5375, length = 604, stop = 7.5375)
    y = range(start = -13.16524, length = 604, stop = 13.16524)
    z = range(start = -0.5, length = 2, stop = 0.5)

    signal = Hist3(x, y, z)
    h = Hist3(x, y, z)

    hX = collect(x)
    kX = collect(y)
    lX = collect(z)
    doctest = MiniVATES.MDNorm(hX, kX, lX)
    maxIx = MiniVATES.maxIntersections(doctest)

    # rot file
    rotFile = benzil_event_nxs_prefix * "0_extra_params.hdf5"
    exData = MiniVATES.loadExtrasData(rotFile)
    m_W = SquareMatrix3c([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])
    transforms2 = MiniVATES.makeTransforms(exData, m_W)

    # sa file
    saData = MiniVATES.loadSolidAngleData(benzil_sa_nxs_file)

    haveSA = true

    # flux file
    fluxData = MiniVATES.loadFluxData(benzil_flux_nxs_file)
    ndets = fluxData.ndets

    # event file
    eventFile = benzil_event_nxs_prefix * "0_BEFORE_MDNorm.nxs"
    eventData = MiniVATES.loadEventData(eventFile)

    intersections = [PreallocVector(Vector{Crd4}(undef, maxIx)) for i = 1:ndets]
    iPerm = [PreallocVector([n for n = 1:maxIx]) for i = 1:ndets]
    yValues = [PreallocVector(Vector{ScalarType}(undef, maxIx)) for i = 1:ndets]

    println("file numbers: (", benzil_event_nxs_min, ", ", benzil_event_nxs_max, ")")
    for file_num = benzil_event_nxs_min:benzil_event_nxs_max
        println("file ", file_num, " / ", )
        fNumStr = string(file_num)
        rotFile = benzil_event_nxs_prefix * fNumStr * "_extra_params.hdf5"
        let extrasWS = MiniVATES.ExtrasWorkspace(rotFile)
            exData.rotMatrix = MiniVATES.getRotationMatrix(extrasWS)
        end

        eventFile = benzil_event_nxs_prefix * fNumStr * "_BEFORE_MDNorm.nxs"
        let eventWS = MiniVATES.EventWorkspace(eventFile)
            eventData.protonCharge = MiniVATES.getProtonCharge(eventWS)
            @time eventData.eventsCtnr, eventData.events =
                MiniVATES.updateEvents!(eventData.eventsCtnr, eventWS)
        end

        transforms = MiniVATES.makeRotationTransforms(exData, m_W)

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

        @time MiniVATES.binEvents!(h, eventData.events, transforms2)
    end

    open("meow.txt", "w") do fio
        outWts = MiniVATES.binweights(h)
        meowWts = MiniVATES.binweights(signal)
        nbins = MiniVATES.nbins(signal)
        for j = 1:nbins[2]
            for i = 1:nbins[1]
                println(fio, outWts[i, j, 1] / meowWts[i, j, 1])
            end
        end
    end
end
