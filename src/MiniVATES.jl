__precompile__(false)
module MiniVATES

import JACC
import Pkg

@static if endswith(JACC.JACCPreferences.backend, "cuda")
    if !haskey(Pkg.project().dependencies, "CUDA")
        # @TODO Julia Pkg.add will add target = :weakdeps in later versions
        Pkg.add("CUDA")
        import CUDA
        CUDA.set_runtime_version!(local_toolkit=true)
    end
    import CUDA
    import CUDA.unsafe_free!
    println("Using CUDA backend for JACC")
elseif endswith(JACC.JACCPreferences.backend, "amdgpu")
    if !haskey(Pkg.project().dependencies, "AMDGPU")
        Pkg.add(; name = "AMDGPU", version = "v0.8.11")
        import AMDGPU
    end
    import AMDGPU
    import AMDGPU.unsafe_free!
    println("Using AMDGPU backend for JACC")
else
    function unsafe_free!(arr) end
    println("Using Threads backend for JACC")
end

import MPI

function __init__()
    # - Initialize here instead of main so that the MPI context can be available
    # for tests.
    # - Conditional allows for case when external users (or tests) have already
    # initialized an MPI context.
    if !MPI.Initialized()
        MPI.Init()
    end
end

include("Util.jl")
include("PreallocArrays.jl")
include("Sort.jl")
include("Hist.jl")
include("BinMD.jl")
include("Load.jl")
include("MDNorm.jl")
include("BinSeries.jl")

end
