include("test_data_constants.jl")
include("common.jl")

import MiniVATES
import JACC
import Atomix
import Adapt
import Test: @testset, @test

@testset "BinMD" begin

    # M = 24
    # N = 291971545
    # tv = Adapt.adapt_structure(JACC.Array, zeros(Int, M))
    # JACC.parallel_for((M, N), (n, i, tv) -> begin
    #     Atomix.@atomic tv[n] += 1
    # end, tv)
    # tvh = Array(tv)
    # for n = 1:M
    #     @test tvh[n] == N
    # end

    x = range(start = -7.5375, length = 604, stop = 7.5375)
    y = range(start = -13.16524, length = 604, stop = 13.16524)
    z = range(start = -0.5, length = 2, stop = 0.5)
    h = MiniVATES.Hist3(x, y, z)
    hf = MiniVATES.Hist3(x, y, z)

    rotFile = benzil_event_nxs_prefix * "0_extra_params.hdf5"
    eventFile = benzil_event_nxs_prefix * "0_BEFORE_MDNorm.nxs"
    fastEventFile = benzil_event_nxs_prefix * "0_events.nxs"
    if length(ARGS) > 0
        base_name = last(ARGS)
        rotFile = base_name * "_extra_params.hdf5"
        eventFile = base_name * "_BEFORE_MDNorm.nxs"
        fastEventFile = base_name * "_events.nxs"
    end

    transforms2 = MiniVATES.makeTransforms(MiniVATES.loadExtrasData(rotFile))
    events = MiniVATES.getEvents(MiniVATES.EventWorkspace(eventFile))
    fastEvents = MiniVATES.getEvents(MiniVATES.FastEventWorkspace(fastEventFile))
    fastWeights = MiniVATES.getWeights(MiniVATES.FastEventWorkspace(fastEventFile))

    # try
        @time MiniVATES.binEvents!(h, events, transforms2)
        @time MiniVATES.binEvents!(hf, fastEvents, fastWeights, transforms2)
    # catch err
    #     code_warntype(err; interactive = true)
    # end
    MiniVATES.reset!(h)
    MiniVATES.reset!(hf)
    @time MiniVATES.binEvents!(h, events, transforms2)
    @time MiniVATES.binEvents!(hf, fastEvents, fastWeights, transforms2)
    @test binweights(h) == binweights(hf)
    write_cat(h)
end
