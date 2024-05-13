module MiniVATES

import JACC
import Pkg

@static if endswith(JACC.JACCPreferences.backend, "cuda")
    # @TODO Julia Pkg.add will add target = :weakdeps in later versions
    Pkg.add("CUDA")
    import CUDA
    println("Using CUDA as back end")
    CUDA.set_runtime_version!(local_toolkit=true)
elseif endswith(JACC.JACCPreferences.backend, "amdgpu")
    Pkg.add("AMDGPU")
    import AMDGPU
    println("Using AMDGPU as back end")
end

include("Util.jl")
include("Hist.jl")
include("MDNorm.jl")

end 
