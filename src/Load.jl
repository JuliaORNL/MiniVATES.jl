import HDF5

struct ExtrasWorkspace
    file::HDF5.File

    ExtrasWorkspace(filename::AbstractString) = new(HDF5.h5open(filename, "r"))
end

@inline function getSkipDets(ws::ExtrasWorkspace)
    return read(ws.file["skip_dets"])
end

@inline function getRotationMatrix(ws::ExtrasWorkspace)
    return transpose(SquareMatrix3c(read(ws.file["expinfo_0"]["goniometer_0"])))
end

@inline function getSymmMatrices(ws::ExtrasWorkspace)
    symm = Vector{SquareMatrix3c}()
    symmGroup = ws.file["symmetryOps"]
    for i = 1:length(symmGroup)
        push!(symm, transpose(read(symmGroup["op_" * string(i - 1)])))
    end
    return symm
end

@inline function getUBMatrix(ws::ExtrasWorkspace)
    return transpose(SquareMatrix3c(read(ws.file["ubmatrix"])))
end

mutable struct ExtrasData
    skip_dets::Vector{Bool}
    rotMatrix::SquareMatrix3c
    symm::Vector{SquareMatrix3c}
    m_UB::SquareMatrix3c
end

function loadExtrasData(rot_nxs_file::AbstractString)
    let ws = ExtrasWorkspace(rot_nxs_file)
        skip_dets = getSkipDets(ws)
        rotMatrix = getRotationMatrix(ws)
        symm = getSymmMatrices(ws)
        m_UB = getUBMatrix(ws)
        return ExtrasData(skip_dets, rotMatrix, symm, m_UB)
    end
end

@inline function makeRotationTransforms(d::ExtrasData, m_W::SquareMatrix3c)
    return Array1(map(op -> inv(d.rotMatrix * d.m_UB * op * m_W), d.symm))
end

@inline function makeTransforms(d::ExtrasData, m_W::SquareMatrix3c)
    return Array1(map(op -> inv(d.m_UB * op * m_W), d.symm))
end

struct SolidAngleData
    solidAngleWS::Vector{Vector{ScalarType}}
    solidAngDetToIdx::Dict{Int32,SizeType}
end

function loadSolidAngleData(sa_nxs_file::AbstractString)
    HDF5.h5open(sa_nxs_file, "r") do file
        saGroup = file["mantid_workspace_1"]
        saData = read(saGroup["workspace"]["values"])
        dims = size(saData)
        @assert length(dims) == 2
        @assert dims[1] == 1
        dataSize = dims[2]
        readData = saData[1, :]
        solidAngleWS = map(x -> [x], readData)

        dcData = read(saGroup["instrument"]["detector"]["detector_count"])
        dims = size(dcData)
        @assert length(dims) == 1
        @assert dims[1] == dataSize

        solidAngDetToIdx = Dict{Int32,SizeType}()
        detector::Int32 = 1
        idx::SizeType = 1
        for value in dcData
            for i = 1:value
                push!(solidAngDetToIdx, detector => idx)
                detector += 1
            end
            idx += 1
        end

        return SolidAngleData(solidAngleWS, solidAngDetToIdx)
    end
end

struct FluxData
    integrFlux_x::AbstractRange{ScalarType}
    integrFlux_y::Vector{Vector{ScalarType}}
    fluxDetToIdx::Dict{Int32,SizeType}
    ndets::SizeType
end

function loadFluxData(flux_nxs_file::AbstractString)
    HDF5.h5open(flux_nxs_file, "r") do file
        group = file["mantid_workspace_1"]
        readDataX::Vector{ScalarType} = read(group["workspace"]["axis1"])
        dims = size(readDataX)
        @assert length(dims) == 1
        dataSize = dims[1]
        integrFlux_x =
            range(length = dataSize, start = first(readDataX), stop = last(readDataX))
        @assert isapprox(
            integrFlux_x[2] - integrFlux_x[1],
            readDataX[2] - readDataX[1],
            rtol = 1.0e-8,
        )

        readDataY::Matrix{ScalarType} = read(group["workspace"]["values"])
        dims = size(readDataY)
        @assert length(dims) == 2
        @assert dims[2] == 1
        integrFlux_y = [readDataY[:, 1]]

        dcData = read(group["instrument"]["detector"]["detector_count"])
        dims = size(dcData)
        @assert length(dims) == 1
        @assert dims[1] == 1
        fluxDetToIdx = Dict{Int32,SizeType}()
        detector = 1
        idx = 1
        for value in dcData
            for i = 1:value
                push!(fluxDetToIdx, detector => idx)
                detector += 1
            end
            idx += 1
        end

        detGroup = group["instrument"]["physical_detectors"]
        readData = read(detGroup["number_of_detectors"])
        dims = size(readData)
        @assert dims == (1,)
        ndets::SizeType = readData[1]

        return FluxData(integrFlux_x, integrFlux_y, fluxDetToIdx, ndets)
    end
end

struct EventWorkspace
    file::HDF5.File

    EventWorkspace(filename::AbstractString) = new(HDF5.h5open(filename, "r"))
end

@inline function getProtonCharge(ws::EventWorkspace)
    expGroup = ws.file["MDEventWorkspace"]["experiment0"]
    pcData = read(expGroup["logs"]["gd_prtn_chrg"]["value"])
    @assert size(pcData) == (1,)
    protonCharge = pcData[1]
end

@inline function _getEventsDataset(ws::EventWorkspace)
    return ws.file["MDEventWorkspace"]["event_data"]["event_data"]
end

@inline function updateEvents!(events::Matrix{CoordType}, ws::EventWorkspace)
    # FIXME: need to resize based on dataset
    # need events to be Vector{SVector{8}}
    ds = _getEventsDataset(ws)
    dims = HDF5.get_extent_dims(ds)
    @show dims
    # copyto!(events, _getEventsDataset(ws))
    events = getEvents(ws)
end

@inline function getEvents(ws::EventWorkspace)
    return read(_getEventsDataset(ws))
end

mutable struct EventData
    lowValues::Vector{CoordType}
    highValues::Vector{CoordType}
    protonCharge::ScalarType
    thetaValues::Array1{CoordType}
    phiValues::Array1{CoordType}
    detIDs::Vector{SizeType}
    events::Matrix{CoordType}
end

function loadEventData(event_nxs_file::AbstractString)
    let ws = EventWorkspace(event_nxs_file)
        expGroup = ws.file["MDEventWorkspace"]["experiment0"]
        lowValues = read(expGroup["logs"]["MDNorm_low"]["value"])
        dims = size(lowValues)
        @assert length(dims) == 1
        dataSize = dims[1]
        highValues = read(expGroup["logs"]["MDNorm_high"]["value"])
        dims = size(highValues)
        @assert length(dims) == 1
        @assert dims[1] == dataSize

        protonCharge = getProtonCharge(ws)

        detGroup = expGroup["instrument"]["physical_detectors"]
        thetaData::Array1{CoordType} = read(detGroup["polar_angle"])
        dims = size(thetaData)
        @assert length(dims) == 1
        @assert dims[1] == dataSize

        phiData::Array1{CoordType} = read(detGroup["azimuthal_angle"])
        dims = size(phiData)
        @assert length(dims) == 1
        @assert dims[1] == dataSize

        thetaValues = thetaData
        phiValues = phiData
        JACC.parallel_for(
            dataSize,
            (i, thetaValues, phiValues) -> begin
                thetaValues[i] = deg2rad(thetaValues[i])
                phiValues[i] = deg2rad(phiValues[i])
            end,
            thetaValues,
            phiValues,
        )

        detIDs = read(detGroup["detector_number"])

        events = getEvents(ws)

        return EventData(
            lowValues,
            highValues,
            protonCharge,
            thetaValues,
            phiValues,
            detIDs,
            events,
        )
    end
end
