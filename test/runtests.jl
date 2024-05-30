function create_link(linkPath::AbstractString)
    if length(ARGS) != 1
        error("Please provide directory path for large test data files")
    end
    lfpath = ARGS[1]
    symlink(joinpath(lfpath, basename(linkPath)), linkPath)
end

const test_data_dir = joinpath(@__DIR__, "data")

const calc_intersections_file = joinpath(test_data_dir, "corelli_demo.txt")
const rot_nxs_file = joinpath(test_data_dir, "CORELLI_29782_rotations.hdf5")
const use_dets_file = joinpath(test_data_dir, "use_det.txt")

const sa_nxs_file = joinpath(test_data_dir, "SingleCrystalDiffuseReduction_SA.nxs")
if !isfile(sa_nxs_file)
    create_link(sa_nxs_file)
end

const flux_nxs_file = joinpath(test_data_dir, "SingleCrystalDiffuseReduction_Flux.nxs")
if !isfile(flux_nxs_file)
    create_link(flux_nxs_file)
end

const event_nxs_file = joinpath(test_data_dir, "CORELLI_29782_Before_MDNorm.nxs")
if !isfile(event_nxs_file)
    create_link(event_nxs_file)
end

const norm_nxs_file = joinpath(test_data_dir, "CORELLI_29782_After_MDNorm_symm.nxs")
if !isfile(norm_nxs_file)
    create_link(norm_nxs_file)
end

include("test_MDNorm.jl")
