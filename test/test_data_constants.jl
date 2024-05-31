
const test_data_dir = joinpath(@__DIR__, "data")

const calc_intersections_file = joinpath(test_data_dir, "corelli_demo.txt")
const rot_nxs_file = joinpath(test_data_dir, "CORELLI_29782_rotations.hdf5")
const use_dets_file = joinpath(test_data_dir, "use_det.txt")

const sa_nxs_file = joinpath(test_data_dir, "SingleCrystalDiffuseReduction_SA.nxs")
const flux_nxs_file = joinpath(test_data_dir, "SingleCrystalDiffuseReduction_Flux.nxs")
const event_nxs_file = joinpath(test_data_dir, "CORELLI_29782_Before_MDNorm.nxs")
const norm_nxs_file = joinpath(test_data_dir, "CORELLI_29782_After_MDNorm_symm.nxs")

const benzil_data_dir = joinpath(test_data_dir, "benzil")
const benzil_event_nxs_prefix = joinpath(benzil_data_dir, "CORELLI_")
const benzil_sa_nxs_file = joinpath(benzil_data_dir, "SolidAngle20160720NoCC.nxs")
const benzil_flux_nxs_file = joinpath(benzil_data_dir, "Spectrum20160720NoCC.nxs")
const benzil_event_nxs_min = 0
const benzil_event_nxs_max = 35

const garnet_data_dir = joinpath(test_data_dir, "garnet")
const garnet_rot_nxs_file = joinpath(garnet_data_dir, "CORELLI_extra_params.hdf5")
const garnet_flux_nxs_file = joinpath(garnet_data_dir, "flux_2p5-8.nxs")
const garnet_sa_nxs_file = joinpath(garnet_data_dir, "solid_angle_2p5-8.nxs")
const garnet_event_nxs_file = joinpath(garnet_data_dir, "CORELLI_BEFORE_MDNorm.nxs")
