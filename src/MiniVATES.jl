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
    println("Using CUDA backend for JACC")
elseif endswith(JACC.JACCPreferences.backend, "amdgpu")
    Pkg.add("AMDGPU")
    import AMDGPU
    println("Using AMDGPU backend for JACC")
end

include("Util.jl")
include("Sort.jl")
include("Hist.jl")
include("BinMD.jl")
include("Load.jl")
include("MDNorm.jl")

end
