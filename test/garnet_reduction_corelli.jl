include("test_data_constants.jl")
include("common.jl")

import MiniVATES
import MiniVATES: binweights, reset!
import HDF5
import Test: @testset, @test

@testset "garnet_reduction_corelli" begin
    x = range(start = -10.05, length = 202, stop = 10.05)
    y = range(start = -10.05, length = 202, stop = 10.05)
    z = range(start = -10.05, length = 202, stop = 10.05)

    signal = MiniVATES.Hist3(x, y, z)
    h = MiniVATES.Hist3(x, y, z)

    # rot file
    @show garnet_rot_nxs_file
    exData = MiniVATES.loadExtrasData(garnet_rot_nxs_file)

    # sa file
    saData = MiniVATES.loadSolidAngleData(garnet_sa_nxs_file)

    # flux file
    fluxData = MiniVATES.loadFluxData(garnet_flux_nxs_file)

    # event file
    eventData = MiniVATES.loadEventData(garnet_event_nxs_file)

    transforms = MiniVATES.makeRotationTransforms(exData)
    transforms2 = MiniVATES.makeTransforms(exData)

    doctest = MiniVATES.MDNorm(signal, exData)

    @time doctest(saData, fluxData, eventData, signal, transforms)
    reset!(signal)
    @time doctest(saData, fluxData, eventData, signal, transforms)

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

    @time MiniVATES.binEvents!(h, eventData.events, transforms2)
    reset!(h)
    @time MiniVATES.binEvents!(h, eventData.events, transforms2)

    write_cat(dims, signal, h)
end
