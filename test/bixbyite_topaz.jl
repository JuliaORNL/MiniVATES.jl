include("test_data_constants.jl")
include("common.jl")

import MiniVATES
import MiniVATES: Hist3, C3

import MPI

dur = @elapsed begin
    x = range(start = -16.0, length = 602, stop = 16.0)
    y = range(start = -16.0, length = 602, stop = 16.0)
    # z = range(start = -16.0, length = 602, stop = 16.0)
    z = range(start = -0.1, length = 2, stop = 0.1)

    extras_events_files = Vector{NTuple{2,AbstractString}}()
    for file_num = bixbyite_event_nxs_min:bixbyite_event_nxs_max
        fNumStr = string(file_num)
        exFile = bixbyite_event_nxs_prefix * fNumStr * "_extra_params.hdf5"
        eventFile = bixbyite_event_nxs_prefix * fNumStr * "_BEFORE_MDNorm.nxs"
        push!(extras_events_files, (exFile, eventFile))
    end

    # MiniVATES.verbose()
    signal, h = MiniVATES.binSeries!(
        MiniVATES.parse_args(ARGS),
        # MiniVATES.parse_args(["--partition=histogram"])
        (x,y,z),
	bixbyite_sa_nxs_file,
        bixbyite_flux_nxs_file,
        extras_events_files,
        C3[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    )

end

write_rank_cat(signal, h)

