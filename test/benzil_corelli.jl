include("test_data_constants.jl")

import MiniVATES
import MiniVATES: ScalarType, CoordType
import MiniVATES: Hist3, SquareMatrix3c, Crd4, PreallocVector
import Test: @testset

@testset "benzil_corelli" begin
    x = range(start = -7.5375, length = 604, stop = 7.5375)
    y = range(start = -13.16524, length = 604, stop = 13.16524)
    z = range(start = -0.5, length = 2, stop = 0.5)

    signal = Hist3(x, y, z)
    h = Hist3(x, y, z)

    # rot file
    rotFile = benzil_event_nxs_prefix * "0_extra_params.hdf5"
    exData = MiniVATES.loadExtrasData(rotFile)

    # sa file
    saData = MiniVATES.loadSolidAngleData(benzil_sa_nxs_file)

    haveSA = true

    # flux file
    fluxData = MiniVATES.loadFluxData(benzil_flux_nxs_file)

    # event file
    eventFile = benzil_event_nxs_prefix * "0_BEFORE_MDNorm.nxs"
    eventData = MiniVATES.loadEventData(eventFile)

    m_W = SquareMatrix3c([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])
    transforms2 = MiniVATES.makeTransforms(exData, m_W)

    doctest = MiniVATES.MDNorm(signal, exData)

    println("file numbers: ", benzil_event_nxs_min, "-", benzil_event_nxs_max)
    for file_num = benzil_event_nxs_min:benzil_event_nxs_max
        println("file #", file_num)
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

        @time doctest(saData, fluxData, eventData, signal, transforms)

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
