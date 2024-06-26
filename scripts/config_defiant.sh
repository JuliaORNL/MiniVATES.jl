#!/bin/bash

# Change the first 3 lines appropriately
PROJ_DIR=/lustre/polis/csc266/scratch/4pf
export JULIA_DEPOT_PATH=$PROJ_DIR/julia_depot
MV_DIR=$PROJ_DIR/GrayScott.jl

# remove existing generated Manifest.toml
rm -f $MV_DIR/Manifest.toml
rm -f $MV_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# load required modules
module load PrgEnv-cray-amd
module load cray-mpich
module load julia

# Required to point at underlying modules above
export JULIA_AMDGPU_DISABLE_ARTIFACTS=1

# MPIPreferences to use Cray's MPICH
julia --project=$MV_DIR -e 'using MPIPreferences; MPIPreferences.use_system_binary(; library_names=["libmpi_cray"], mpiexec="srun")'

# Regression being fixed with CUDA v4.0.0. CUDA.jl does lazy loading for portability to systems without NVIDIA GPUs
# julia --project=$MV_DIR -e 'using Pkg; Pkg.add(name="CUDA", version="v3.13.1")' 
# Adds a custom branch in case the development version is needed (for devs to test new features)
julia --project=$MV_DIR -e 'using Pkg; Pkg.add("AMDGPU")'
julia --project=$MV_DIR -e 'using JACC; JACC.JACCPreferences.set_backend("amdgpu")'
julia --project -e 'using JACC.JACCPreferences; JACCPreferences.set_backend("CUDA")'

# Instantiate the project by installing packages in Project.toml
julia --project=$MV_DIR -e 'using Pkg; Pkg.instantiate()'

# Verify the packages are installed correctly
julia --project=$MV_DIR -e 'using Pkg; Pkg.build()'
julia --project=$MV_DIR -e 'using Pkg; Pkg.precompile()'
