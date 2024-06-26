import MPI
using Printf

@inline function getRankRange(N::Integer)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    size = MPI.Comm_size(comm)
    count = trunc(Int, N / size)
    remainder = trunc(Int, N % size)
    if rank < remainder
        # The first 'remainder' ranks get 'count + 1' tasks each
        start = rank * (count + 1)
        stop = start + count
    else
        # The remaining 'size - remainder' ranks get 'count' task each
        start = rank * count + remainder
        stop = start + (count - 1)
    end

    return (start + 1, stop + 1)
end

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

    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    commSize = MPI.Comm_size(MPI.COMM_WORLD)
    start, stop = getRankRange(length(eventFilePairs))

    if MiniVATES.be_verbose
        if rank == 0
            println("number of files: ", length(eventFilePairs))
        end
        @show rank, start, stop
    end

    for fi = start:stop
        exFile, eventFile = eventFilePairs[fi]
        let extrasWS = ExtrasWorkspace(exFile)
            exData.rotMatrix = getRotationMatrix(extrasWS)
        end

        updateEventsTime = nothing
        let eventWS = EventWorkspace(eventFile)
            eventData.protonCharge = getProtonCharge(eventWS)
            dur = @timed begin
                updateEvents!(eventData, eventWS)
            end
            updateEventsTime = dur.time
        end

        transforms = makeRotationTransforms(exData)

        dur = @timed begin
            mdNorm!(signal, mdn, saData, fluxData, eventData, transforms)
        end
        mdNormTime = dur.time

        dur = @timed begin
            binEvents!(eventsHist, eventData.events, transforms2)
        end
        binEventsTime = dur.time

        for r = 0:(commSize-1)
            if rank == r
                println(
                    "rank: ",
                    lpad(rank, 2),
                    "; fi: ",
                    lpad(fi, 3),
                    "; updateEvents: ",
                    @sprintf("%3.6f", updateEventsTime),
                    ", mdNorm: ",
                    @sprintf("%3.6f", mdNormTime),
                    ", binEvents: ",
                    @sprintf("%3.6f", binEventsTime),
                )
            end
            MPI.Barrier(MPI.COMM_WORLD)
        end
    end
    # MPI.Barrier(MPI.COMM_WORLD)

    return (signal, eventsHist)
end
