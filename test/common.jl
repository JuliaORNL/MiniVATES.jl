import MiniVATES: Hist3, binweights, nbins, SignalType

function write_cat(signal::Union{Hist3, Nothing}, h::Hist3)
    dims = nbins(h)
    open("purr.txt", "w") do fio
        println(fio, dims[1])
        println(fio, dims[2])
    end

    outWts = Array(MiniVATES.binweights(h))
    meowWts = ones(SignalType, dims)
    if signal isa Hist3
        meowWts = Array(MiniVATES.binweights(signal))
    end

    open("meow.txt", "w") do fio
        for j = 1:dims[2]
            for i = 1:dims[1]
                println(fio, outWts[i, j, 1] / meowWts[i, j, 1])
            end
        end
    end
end

write_cat(h::Hist3) = write_cat(nothing, h)