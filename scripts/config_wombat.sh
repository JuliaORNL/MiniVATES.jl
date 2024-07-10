#!/bin/bash

MV_DIR=/ccsopen/home/4pf/MiniVATES.jl

module load nvhpc
# module load hpcx/hpcx-stack
# module load hpcx/hpcx-ompi

export CUDA_HOME=/sw/wombat/Nvidia_HPC_SDK/Linux_aarch64/24.5/cuda/12.4

export PATH=$PATH:/ccsopen/home/4pf/.julia/bin

# rm -f $MV_DIR/Manifest.toml
# rm -f $MV_DIR/LocalPreferences.toml

julia --project=$MV_DIR -e 'using Pkg; Pkg.instantiate()'

# MPI + CUDA
julia --project=$MV_DIR -e ' \
    using Pkg; \
    if !haskey(Pkg.project().dependencies, "MPIPreferences"); \
        Pkg.add("MPIPreferences"); \
    end; \
    using MPIPreferences; \
    MPIPreferences.use_jll_binary(); \
    using MPI; \
    MPI.install_mpiexecjl(force = true); \
    if !haskey(Pkg.project().dependencies, "CUDA"); \
        Pkg.add("CUDA"); \
    end; \
    using CUDA; \
    CUDA.set_runtime_version!(v"12.4"; local_toolkit=true); \
    '

# JACC
julia --project=$MV_DIR -e ' \
    using Pkg; \
    jaccInfo = Pkg.dependencies()[Pkg.project().dependencies["JACC"]]; \
    @show jaccInfo.git_revision; \
    if jaccInfo.git_revision != "main"; \
        Pkg.add(; name="JACC", rev = "main"); \
    end; \
    using JACC; \
    JACC.JACCPreferences.set_backend("cuda"); \
    '

# Verify the packages are installed correctly
julia --project=$MV_DIR -e 'using Pkg; Pkg.instantiate()'
julia --project=$MV_DIR -e ' \
    using Pkg; \
    Pkg.status(); \
    using CUDA; \
    CUDA.versioninfo(); \
    '
