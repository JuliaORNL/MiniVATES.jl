function create_link(linkPath::AbstractString)
    if length(ARGS) != 1
        error("Please provide directory path for large test data files")
    end
    lfpath = ARGS[1]
    symlink(joinpath(lfpath, basename(linkPath)), linkPath)
end

include("test_data_constants.jl")

if !isfile(sa_nxs_file)
    create_link(sa_nxs_file)
end

if !isfile(flux_nxs_file)
    create_link(flux_nxs_file)
end

if !isfile(event_nxs_file)
    create_link(event_nxs_file)
end

if !isfile(norm_nxs_file)
    create_link(norm_nxs_file)
end

if !ispath(benzil_data_dir)
    create_link(benzil_data_dir)
end

if !ispath(garnet_data_dir)
    create_link(garnet_data_dir)
end
