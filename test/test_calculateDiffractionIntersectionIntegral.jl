include("common.jl")
include("test_data_constants.jl")

import MiniVATES
import MiniVATES: ScalarType, CoordType
import MiniVATES: SquareMatrix3c, Crd4
import MiniVATES: Hist3, binweights, reset!
import Test: @test, @testset
import HDF5

using JACC
using Cthulhu
import Profile
# import StatProfilerHTML: statprofilehtml

@testset "calculateDiffractionIntersectionIntegral" begin
    x = range(-10.0, length = 201, stop = 10.0)
    y = range(-10.0, length = 201, stop = 10.0)
    z = range(-0.1, length = 2, stop = 0.1)

    signal = Hist3(x, y, z)

    # extras file
    exData = MiniVATES.loadExtrasData(rot_nxs_file)

    # sa file
    saData = MiniVATES.loadSolidAngleData(sa_nxs_file)

    # flux file
    fluxData = MiniVATES.loadFluxData(flux_nxs_file)

    # event file
    eventData = MiniVATES.loadEventData(event_nxs_file)

    MiniVATES.set_m_W!(exData, [1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0 0.0 1.0])
    transforms = MiniVATES.makeRotationTransforms(exData)
    transforms2 = MiniVATES.makeTransforms(exData)

    doctest = MiniVATES.MDNorm(x, y, z, exData)

    try
        @time MiniVATES.mdNorm!(signal, doctest, saData, fluxData, eventData, transforms)
    catch err
        code_warntype(err; interactive = true)
    end
    reset!(signal)
    @time MiniVATES.mdNorm!(signal, doctest, saData, fluxData, eventData, transforms)
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
    ref_max = 0.0
    signalWts = Core.Array(binweights(signal))
    max_signal = maximum(signalWts)
    for i = 1:dims[2]
        for j = 1:dims[1]
            @test isapprox(data2d[i, j], signalWts[i, j, 1], atol = 1.0e+6)
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

    write_cat(signal, h)
end
