import HDF5

struct RotationData
    rotMatrix::SquareMatrix3c
    symm::Vector{SquareMatrix3c}
    m_UB::SquareMatrix3c
    m_W::SquareMatrix3c
end

function loadRotationData(rot_nxs_file::AbstractString)
    HDF5.h5open(rot_nxs_file, "r") do file
        rotMatrix = transpose(SquareMatrix3c(read(file["expinfo_0"]["goniometer_0"])))

        symm = Vector{SquareMatrix3c}()
        symmGroup = file["symmetryOps"]
        for i = 1:length(symmGroup)
            push!(symm, transpose(read(symmGroup["op_" * string(i - 1)])))
        end

        m_UB = transpose(SquareMatrix3c(read(file["ubmatrix"])))
        m_W = SquareMatrix3c([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])

        return RotationData(rotMatrix, symm, m_UB, m_W)
    end
end

@inline function makeTransforms(rd::RotationData)
    return Array1(map(op -> inv(rd.rotMatrix * rd.m_UB * op * rd.m_W), rd.symm))
end

struct SAData
    solidAngleWS::Vector{Vector{ScalarType}}
    solidAngDetToIdx::Dict{Int32,SizeType}
end

function loadSAData(sa_nxs_file::AbstractString)
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

        return SAData(solidAngleWS, solidAngDetToIdx)
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

mutable struct EventsData
    lowValues::Vector{CoordType}
    highValues::Vector{CoordType}
    protonCharge::ScalarType
    thetaValues::Vector{CoordType}
    phiValues::Vector{CoordType}
    detIDs::Vector{SizeType}
    events::Matrix{CoordType}
end

function loadEventsData(event_nxs_file::AbstractString)
    HDF5.h5open(event_nxs_file, "r") do file
        expGroup = file["MDEventWorkspace"]["experiment0"]
        lowValues = read(expGroup["logs"]["MDNorm_low"]["value"])
        dims = size(lowValues)
        @assert length(dims) == 1
        dataSize = dims[1]
        highValues = read(expGroup["logs"]["MDNorm_high"]["value"])
        dims = size(highValues)
        @assert length(dims) == 1
        @assert dims[1] == dataSize

        pcData = read(expGroup["logs"]["gd_prtn_chrg"]["value"])
        @assert size(pcData) == (1,)
        protonCharge = pcData[1]

        detGroup = expGroup["instrument"]["physical_detectors"]
        thetaData::Vector{CoordType} = read(detGroup["polar_angle"])
        dims = size(thetaData)
        @assert length(dims) == 1
        @assert dims[1] == dataSize
        thetaValues = map(deg2rad, thetaData)

        phiData::Vector{CoordType} = read(detGroup["azimuthal_angle"])
        dims = size(phiData)
        @assert length(dims) == 1
        @assert dims[1] == dataSize
        phiValues = map(deg2rad, phiData)

        detIDs = read(detGroup["detector_number"])

        evGroup = file["MDEventWorkspace"]["event_data"]
        events = read(evGroup["event_data"])

        return EventsData(
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
