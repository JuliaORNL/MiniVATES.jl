import MiniVATES
import MiniVATES: SizeType, ScalarType, CoordType, SquareMatrix3, V3, Vec4
import MiniVATES: Hist3, atomic_push!, binweights
import Test: @test, @testset
import HDF5

using Profile
using StatProfilerHTML

@testset "calculateIntersections" begin
    x = range(-10.0, length = 201, stop = 10.0)
    y = range(-10.0, length = 201, stop = 10.0)
    z = range(-0.1, length = 2, stop = 0.1)

    histogram = Hist3(x, y, z)

    hX = collect(x)
    kX = collect(y)
    lX = collect(z)

    doctest = MiniVATES.MDNorm(hX, kX, lX)

    open(calc_intersections_file) do f
        line = split(strip(readline(f)), ' ')
        transform = transpose(SquareMatrix3(parse.(Float64, line)))

        ndets = parse(Int, readline(f))

        intersections = Vector{Vec4}()
        for i = 1:ndets
            line = split(strip(readline(f)), ' ')
            i_f = parse(Int, line[1])
            theta = parse(Float64, line[2])
            phi = parse(Float64, line[3])
            lowvalue = parse(Float64, line[4])
            highvalue = parse(Float64, line[5])
            if i != i_f
                i = i_f
            end
            num_intersections = parse(Int, readline(f))
            values = Vector{Vector{Float64}}()
            for i = 1:num_intersections
                line = split(strip(readline(f)), ' ')
                push!(values, parse.(Float64, line))
            end

            MiniVATES.calculateIntersections!(
                doctest,
                histogram,
                theta,
                phi,
                transform,
                lowvalue,
                highvalue,
                intersections,
            )
            @test length(intersections) == num_intersections
            for i = 1:num_intersections
                for j = 1:4
                    @test abs(intersections[i][j] - values[i][j]) < 0.0001
                end
            end
        end
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
    readDataX = read(group["workspace"]["axis1"])
    dims = size(readDataX)
    @test length(dims) == 1
    integrFlux_x =
        range(length = length(readDataX), start = first(readDataX), stop = last(readDataX))
    @test isapprox(
        integrFlux_x[2] - integrFlux_x[1],
        readDataX[2] - readDataX[1],
        rtol = 1.0e-8,
    )

    readDataY = read(group["workspace"]["values"])
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

    readData = read(group["instrument"]["physical_detectors"]["number_of_detectors"])
    dims = size(readData)
    @test length(dims) == 1
    @test dims[1] == 1
    ndets = readData[1]

    close(file)

    # use dets file
    use_dets = parse.(Bool, readlines(use_dets_file))
    for v in use_dets
        @test v == true
    end

    # event file
    file = HDF5.h5open(event_nxs_file, "r")
    expGroup = file["MDEventWorkspace"]["experiment0"]
    lowValues = read(expGroup["logs"]["MDNorm_low"]["value"])
    dims = size(lowValues)
    @test length(dims) == 1
    @test dims[1] == 372736
    highValues = read(expGroup["logs"]["MDNorm_high"]["value"])
    dims = size(highValues)
    @test length(dims) == 1
    @test dims[1] == 372736

    pcData = read(expGroup["logs"]["gd_prtn_chrg"]["value"])
    dims = size(pcData)
    @test length(dims) == 1
    @test dims[1] == 1
    protonCharge = pcData[1]

    thetaData = read(expGroup["instrument"]["physical_detectors"]["polar_angle"])
    dims = size(thetaData)
    @test length(dims) == 1
    @test dims[1] == 372736
    thetaValues = map(deg2rad, thetaData)

    phiData = read(expGroup["instrument"]["physical_detectors"]["azimuthal_angle"])
    dims = size(phiData)
    @test length(dims) == 1
    @test dims[1] == 372736
    phiValues = map(deg2rad, phiData)

    detIDs = read(expGroup["instrument"]["physical_detectors"]["detector_number"])
    dims = size(detIDs)
    @test length(dims) == 1
    @test dims[1] == 372736

    events = read(file["MDEventWorkspace"]["event_data"]["event_data"])
    dims = size(events)
    @test length(dims) == 2
    close(file)

    MDDims = 3
    transforms = map(op -> inv(rotMatrix * m_UB * op * m_W), symm)
    nSymm = length(symm)
    intersections = Vector{Vec4}()

    function kernel1()
        signal = Hist3(x, y, z)
        @time begin
            for n = 1:nSymm
                for i = 1:ndets
                    @inbounds begin
                        if !use_dets[i]
                            continue
                        end

                        detID = detIDs[i]
                        wsIdx = get(fluxDetToIdx, detID, nothing)
                        if wsIdx == nothing
                            continue
                        end

                        MiniVATES.calculateIntersections!(
                            doctest,
                            signal,
                            thetaValues[i],
                            phiValues[i],
                            transforms[n],
                            lowValues[i],
                            highValues[i],
                            intersections,
                        )

                        if isempty(intersections)
                            continue
                        end

                        solid = protonCharge
                        if haveSA
                            solidAngleFactor = solidAngleWS[solidAngDetToIdx[detID]][1]
                            solid *= solidAngleFactor
                        end
                        xValues, yValues =
                            MiniVATES.calculateDiffractionIntersectionIntegral(
                                intersections,
                                integrFlux_x,
                                integrFlux_y,
                                wsIdx,
                            )

                        MiniVATES.calculateSingleDetectorNorm!(
                            doctest,
                            intersections,
                            solid,
                            yValues,
                            signal,
                        )
                    end
                end
            end
        end
    end
    kernel1()
    GC.enable(false)
    kernel1()
    # Profile.clear()
    # @profile kernel1()
    GC.enable(true)

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
    # max_signal = maximum(signal.weights)
    max_signal = 0.0
    ref_max = 0.0
    for i = 1:dims[2]
        for j = 1:dims[1]
            @test isapprox(data2d[i, j], binweights(signal)[i, j, 1], atol = 1.0e+6)
            ref_max = max(ref_max, data2d[j, i])
            max_signal = max(max_signal, binweights(signal)[i, j, 1])
        end
    end
    @test isapprox(max_signal, ref_max, atol = 1.0e+5)

    # h = StatsBase.Histogram((x, y, z), Float64)
    h = Hist3(x, y, z)

    function kernel2()
        @time begin
            @inbounds begin
                for op in transforms
                    for i = 1:size(events)[2]
                        v = op * V3[events[6, i], events[7, i], events[8, i]]
                        atomic_push!(h, v[1], v[2], v[3], events[1, i])
                    end
                end
            end
        end
    end
    kernel2()
    GC.enable(false)
    kernel2()
    # @profile kernel2()
    # statprofilehtml()
    GC.enable(true)

    fio = open("meow.txt", "w")
    for j = 1:dims[2]
        for i = 1:dims[1]
            println(fio, binweights(h)[i, j, 1] / binweights(signal)[i, j, 1])
        end
    end
end
