import MiniVATES
import MiniVATES: SizeType, ScalarType, CoordType, SquareMatrix3, V3, Vec4
import MiniVATES: Hist3, atomic_push!, binweights, reset!, maxIntersections
import MiniVATES: PreallocVector
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

#     open(calc_intersections_file) do f
#         line = split(strip(readline(f)), ' ')
#         transform = transpose(SquareMatrix3(parse.(ScalarType, line)))

#         ndets = parse(Int, readline(f))

#         intersections = Vector{Vec4}()
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

#             MiniVATES.calculateIntersections!(
#                 doctest,
#                 histogram,
#                 theta,
#                 phi,
#                 transform,
#                 lowvalue,
#                 highvalue,
#                 intersections,
#             )
#             @test length(intersections) == num_intersections
#             for i = 1:num_intersections
#                 for j = 1:4
#                     @test abs(intersections[i][j] - values[i][j]) < 0.0001
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

# function kernel1_2D(n, i, t)
#     @inbounds begin
#         if !t.use_dets[i]
#             return
#         end

#         detID = t.detIDs[i]
#         wsIdx = get(t.fluxDetToIdx, detID, nothing)
#         if wsIdx == nothing
#             return
#         end

#         intersections = MiniVATES.calculateIntersections(
#             t.doctest,
#             t.signal,
#             t.thetaValues[i],
#             t.phiValues[i],
#             t.transforms[n],
#             t.lowValues[i],
#             t.highValues[i],
#         )

#         if isempty(intersections)
#             return
#         end

#         yValues = MiniVATES.calculateDiffractionIntersectionIntegral(
#             intersections,
#             t.integrFlux_x,
#             t.integrFlux_y[wsIdx],
#         )

#         solid::ScalarType = t.protonCharge
#         if t.haveSA
#             saIdx = t.solidAngDetToIdx[detID]
#             saFactor = t.solidAngleWS[saIdx][1]
#             solid *= saFactor
#         end
#         MiniVATES.calculateSingleDetectorNorm!(
#             t.doctest,
#             intersections,
#             solid,
#             yValues,
#             t.signal,
#         )
#     end
# end

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
    file = HDF5.h5open(rot_nxs_file, "r")
    rotMatrix = transpose(SquareMatrix3(read(file["expinfo_0"]["goniometer_0"])))

    symm = Vector{SquareMatrix3}()
    symmGroup = file["symmetryOps"]
    for i = 1:length(symmGroup)
        push!(symm, transpose(read(symmGroup["op_" * string(i - 1)])))
    end

    m_UB = transpose(SquareMatrix3(read(file["ubmatrix"])))
    close(file)
    m_W = SquareMatrix3([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])

    # sa file
    file = HDF5.h5open(sa_nxs_file, "r")
    saGroup = file["mantid_workspace_1"]
    saData = read(saGroup["workspace"]["values"])
    dims = size(saData)
    @test length(dims) == 2
    @test dims[1] == 1
    @test dims[2] == 372736
    readData = saData[1, :]
    solidAngleWS = map(x -> [x], readData)

    dcData = read(saGroup["instrument"]["detector"]["detector_count"])
    dims = size(dcData)
    @test length(dims) == 1
    @test dims[1] == 372736

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

    close(file)

    haveSA = true

    # flux file
    file = HDF5.h5open(flux_nxs_file, "r")
    group = file["mantid_workspace_1"]
    readDataX::Vector{ScalarType} = read(group["workspace"]["axis1"])
    dims = size(readDataX)
    @test length(dims) == 1
    integrFlux_x =
        range(length = length(readDataX), start = first(readDataX), stop = last(readDataX))
    @test isapprox(
        integrFlux_x[2] - integrFlux_x[1],
        readDataX[2] - readDataX[1],
        rtol = 1.0e-8,
    )

    readDataY::Matrix{ScalarType} = read(group["workspace"]["values"])
    dims = size(readDataY)
    @test length(dims) == 2
    @test dims[2] == 1
    integrFlux_y = [readDataY[:, 1]]

    dcData = read(group["instrument"]["detector"]["detector_count"])
    dims = size(dcData)
    @test length(dims) == 1
    @test dims[1] == 1

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
    @test length(dims) == 1
    @test dims[1] == 1
    ndets::SizeType = readData[1]

    close(file)

    # use dets file
    use_dets = parse.(Bool, readlines(use_dets_file))
    for v in use_dets
        @test v == true
    end

    # event file
    file = HDF5.h5open(event_nxs_file, "r")
    expGroup = file["MDEventWorkspace"]["experiment0"]
    lowValues::Vector{CoordType} = read(expGroup["logs"]["MDNorm_low"]["value"])
    dims = size(lowValues)
    @test length(dims) == 1
    @test dims[1] == 372736
    highValues::Vector{CoordType} = read(expGroup["logs"]["MDNorm_high"]["value"])
    dims = size(highValues)
    @test length(dims) == 1
    @test dims[1] == 372736

    pcData = read(expGroup["logs"]["gd_prtn_chrg"]["value"])
    dims = size(pcData)
    @test length(dims) == 1
    @test dims[1] == 1
    protonCharge::ScalarType = pcData[1]

    detGroup = expGroup["instrument"]["physical_detectors"]
    thetaData::Vector{CoordType} = read(detGroup["polar_angle"])
    dims = size(thetaData)
    @test length(dims) == 1
    @test dims[1] == 372736
    thetaValues = map(deg2rad, thetaData)

    phiData::Vector{CoordType} = read(detGroup["azimuthal_angle"])
    dims = size(phiData)
    @test length(dims) == 1
    @test dims[1] == 372736
    phiValues = map(deg2rad, phiData)

    detIDs = read(detGroup["detector_number"])
    dims = size(detIDs)
    @test length(dims) == 1
    @test dims[1] == 372736

    evGroup = file["MDEventWorkspace"]["event_data"]
    events::Matrix{CoordType} = read(evGroup["event_data"])
    dims = size(events)
    @test length(dims) == 2
    close(file)

    MDDims = 3
    transforms = JACC.Vector(map(op -> inv(rotMatrix * m_UB * op * m_W), symm))
    nSymm::SizeType = length(symm)
    maxIx = maxIntersections(doctest)
    intersections = [PreallocVector(Vector{Vec4}(undef, maxIx)) for i = 1:ndets]
    yValues = [PreallocVector(Vector{ScalarType}(undef, maxIx)) for i = 1:ndets]
    iPerm = [PreallocVector([n for n = 1:maxIx]) for i = 1:ndets]

    function launch_kernel1()
        @time begin
            for n = 1:nSymm
                JACC.parallel_for(
                    ndets,
                    kernel1_1D,
                    (
                        transform = transforms[n],
                        use_dets = use_dets,
                        intersections,
                        iPerm,
                        yValues,
                        detIDs,
                        fluxDetToIdx,
                        doctest,
                        signal,
                        thetaValues,
                        phiValues,
                        lowValues,
                        highValues,
                        integrFlux_x,
                        integrFlux_y,
                        solidAngDetToIdx,
                        solidAngleWS,
                        protonCharge,
                        haveSA,
                    ),
                )
            end
        end
    end
    launch_kernel1()
    # reset!(signal)
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

    @time begin
        JACC.parallel_for(
            (nSymm, size(events)[2]),
            (n, i, t) -> begin
                op = t.transforms[n]
                v = op * V3[t.events[6, i], t.events[7, i], t.events[8, i]]
                atomic_push!(t.h, v[1], v[2], v[3], t.events[1, i])
            end,
            (h = h, events, transforms),
        )
    end

    fio = open("meow.txt", "w")
    for j = 1:dims[2]
        for i = 1:dims[1]
            println(fio, binweights(h)[i, j, 1] / binweights(signal)[i, j, 1])
        end
    end
end
