import MiniVATES: Hist3, binweights, nbins

function write_cat(signal::Hist3, h::Hist3)
    dims = nbins(signal)

    open("purr.txt", "w") do fio
        println(fio, dims[1])
        println(fio, dims[2])
    end

    outWts = MiniVATES.binweights(h)
    meowWts = MiniVATES.binweights(signal)

    open("meow.txt", "w") do fio
        for j = 1:dims[2]
            for i = 1:dims[1]
                println(fio, outWts[i, j, 1] / meowWts[i, j, 1])
            end
        end
    end
end
