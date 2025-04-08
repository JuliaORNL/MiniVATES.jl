import MiniVATES: Hist3, binweights, nbins, SignalType

function write_cat(signal::Union{Hist3,Nothing}, h::Hist3, file_suffix = "")
    dims = nbins(h)
    open("purr" * file_suffix * ".txt", "w") do fio
        println(fio, dims[1])
        println(fio, dims[2])
    end

    outWts = Core.Array(MiniVATES.binweights(h))
    if signal isa Hist3
        meowWts = Core.Array(MiniVATES.binweights(signal))
    else
        meowWts = ones(SignalType, dims)
    end

    open("meow" * file_suffix * ".txt", "w") do fio
        for i = 1:dims[1]
            for j = 1:dims[2]
                println(fio, outWts[i, j, 1] / meowWts[i, j, 1])
            end
        end
    end
end

write_cat(h::Hist3) = write_cat(nothing, h)

function write_rank_cat(signal::Hist3, h::Hist3)
    rank_str = string(MPI.Comm_rank(MPI.COMM_WORLD))
    write_cat(signal, h, "_" * rank_str)
end
