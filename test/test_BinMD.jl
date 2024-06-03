include("test_data_constants.jl")
include("common.jl")

import MiniVATES
import Test: @testset

@testset "BinMD" begin
    x = range(start = -7.5375, length = 604, stop = 7.5375)
    y = range(start = -13.16524, length = 604, stop = 13.16524)
    z = range(start = -0.5, length = 2, stop = 0.5)
    h = MiniVATES.Hist3(x, y, z)

    rotFile = benzil_event_nxs_prefix * "0_extra_params.hdf5"
    eventFile = benzil_event_nxs_prefix * "0_BEFORE_MDNorm.nxs"
    if length(ARGS) > 0
        base_name = last(ARGS)
        rotFile = base_name * "_extra_params.hdf5"
        eventFile = base_name * "_BEFORE_MDNorm.nxs"
    end

    @show rotFile
    @show eventFile
    transforms2 = MiniVATES.makeTransforms(MiniVATES.loadExtrasData(rotFile))
    events = MiniVATES.getEvents(MiniVATES.EventWorkspace(eventFile))

    # try
        @time MiniVATES.binEvents!(h, events, transforms2)
    # catch err
    #     code_warntype(err; interactive = true)
    # end
    MiniVATES.reset!(h)
    @time MiniVATES.binEvents!(h, events, transforms2)

    write_cat(h)
end
