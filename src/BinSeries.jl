import MPI
using Printf

tmfmt(tm::AbstractFloat) = @sprintf("%3.6fs", tm)

@inline function binSeries!(
    signal::Hist3,
    eventsHist::Hist3,
    mdn::MDNorm,
    saFile::AbstractString,
    fluxFile::AbstractString,
    eventFilePairs::Vector{NTuple{2,AbstractString}},
    m_W::SquareMatrix3c,
)
    saData = loadSolidAngleData(saFile)
    fluxData = loadFluxData(fluxFile)

    exFile, eventFile = first(eventFilePairs)
    exData = loadExtrasData(exFile)
    setExtrasData!(mdn, exData)
    eventData = loadEventData(eventFile)

    set_m_W!(exData, m_W)
    transforms2 = makeTransforms(exData)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    commSize = MPI.Comm_size(comm)
    nFiles = length(eventFilePairs)
    # start, stop = getRankRange(nFiles)
    start, stop = (1, nFiles)

    if MiniVATES.be_verbose
        if rank == 0
            println("number of files: ", nFiles)
        end
        @show rank, start, stop
    end

    updAvgJ = 0.0
    mdnAvgJ = 0.0
    binAvgJ = 0.0

    updAvg = 0.0
    mdnAvg = 0.0
    binAvg = 0.0

    for fi = start:stop
        exFile, eventFile = eventFilePairs[fi]
        let extrasWS = ExtrasWorkspace(exFile)
            exData.rotMatrix = getRotationMatrix(extrasWS)
        end

        updateEventsTime = nothing
        let eventWS = EventWorkspace(eventFile)
            eventData.protonCharge = getProtonCharge(eventWS)
            dur = @elapsed updateEvents!(eventData, eventWS)
            updateEventsTime = dur
            updAvg += dur
        end

        transforms = makeRotationTransforms(exData)

        dur = @elapsed mdNorm!(signal, mdn, saData, fluxData, eventData, transforms)
        mdNormTime = dur
        mdnAvg += dur

        dur = @elapsed binEvents!(eventsHist, eventData.events, transforms2)
        binEventsTime = dur
        binAvg += dur

        for r = 0:(commSize - 1)
            if rank == r
                println(
                    "rank: ",
                    lpad(rank, 2),
                    "; fi: ",
                    lpad(fi, 3),
                    "; updateEvents: ",
                    tmfmt(updateEventsTime),
                    ", mdNorm: ",
                    tmfmt(mdNormTime),
                    ", binEvents: ",
                    tmfmt(binEventsTime),
                )
            end
            MPI.Barrier(comm)
        end

        if fi == start
            updAvgJ = updAvg
            mdnAvgJ = mdnAvg
            binAvgJ = binAvg
            updAvg = 0.0
            mdnAvg = 0.0
            binAvg = 0.0
        end
    end
    sumJ = [updAvgJ, mdnAvgJ, binAvgJ]
    MPI.Reduce!(sumJ, MPI.SUM, 0, comm)
    sum = [updAvg, mdnAvg, binAvg]
    MPI.Reduce!(sum, MPI.SUM, 0, comm)
    if rank == 0
        avgJ = sumJ ./ commSize
        avg = sum ./ nFiles
        println("Averages:")
        println("    updateEvents (JIT): ", tmfmt(avgJ[1]))
        println("    updateEvents:       ", tmfmt(avg[1]))
        println("    mdNorm (JIT):       ", tmfmt(avgJ[2]))
        println("    mdNorm:             ", tmfmt(avg[2]))
        println("    binEvents (JIT):    ", tmfmt(avgJ[3]))
        println("    binEvents:          ", tmfmt(avg[3]))
    end

    return (signal, eventsHist)
end

function mergeHistogramToRootProcess!(hist::Hist3)
    dur = @elapsed begin
        weights = Core.Array(binweights(hist))
        MPI.Reduce!(weights, MPI.SUM, 0, MPI.COMM_WORLD)
        ret = Hist3(
            edges(hist),
            nbins(hist),
            origin(hist),
            boxLength(hist),
            adapt_structure(JACCArray, weights),
        )
    end
    if MPI.Comm_rank(MPI.COMM_WORLD) == 0
        println("Reduce: ", tmfmt(dur))
    end
    return ret
end
