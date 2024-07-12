#!/bin/bash

# Change the first 3 lines appropriately
PROJ_DIR=/lustre/polis/csc266/scratch/4pf
export JULIA_DEPOT_PATH=$PROJ_DIR/julia_depot
MV_DIR=$PROJ_DIR/MiniVATES.jl
echo $MV_DIR

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-cray-amd
module load cray-mpich
module load julia

# remove existing generated Manifest.toml
rm -f $MV_DIR/Manifest.toml
rm -f $MV_DIR/LocalPreferences.toml

# Required to point at underlying modules above
export JULIA_AMDGPU_DISABLE_ARTIFACTS=1

julia --project=$MV_DIR -e 'using Pkg; Pkg.instantiate()'

# MPI + AMDGPU
julia --project=$MV_DIR -e ' \
    using Pkg; \
    if !haskey(Pkg.project().dependencies, "MPIPreferences"); \
        Pkg.add("MPIPreferences"); \
    end; \
    using MPIPreferences; \
    MPIPreferences.use_system_binary(mpiexec="srun", vendor="cray"); \
    if !haskey(Pkg.project().dependencies, "AMDGPU"); \
        Pkg.add(; name="AMDGPU", version = "v0.8.11"); \
    end; \
    '

# JACC
julia --project=$MV_DIR -e ' \
    using Pkg; \
    jaccInfo = Pkg.dependencies()[Pkg.project().dependencies["JACC"]]; \
    if jaccInfo.git_revision != "main"; \
        Pkg.add(; name="JACC", rev = "main"); \
    end; \
    using JACC; \
    JACC.JACCPreferences.set_backend("amdgpu"); \
    '

# Verify the packages are installed correctly
julia --project=$MV_DIR -e 'using Pkg; Pkg.instantiate()'
julia --project=$MV_DIR -e ' \
    using Pkg; \
    Pkg.status(); \
    '
