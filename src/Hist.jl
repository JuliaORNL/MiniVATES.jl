import FHist
import Atomix
import Base: @propagate_inbounds

###

struct Hist3_FHist
    impl::FHist.Hist3D

    function Hist3_FHist(x::AbstractArray, y::AbstractArray, z::AbstractArray)
        new(FHist.Hist3D(counttype = SignalType, binedges = (x, y, z)))
    end
end

@inline function atomic_push!(h::Hist3_FHist, x, y, z, wt)
    FHist.atomic_push!(h.impl, x, y, z, wt)
end

@inline binweights(h::Hist3_FHist) = FHist.bincounts(h.impl)

@inline function binindex1d(r, val)
    if !(first(r) <= val)
        return 0
    elseif !(val <= last(r))
        return length(r)
    else
        return searchsortedlast(r, val)
    end
end

@inline function binindex(h::Hist3_FHist, x, y, z)
    rx, ry, rz = FHist.binedges(h.impl)
    return (binindex1d(rx, x), binindex1d(ry, y), binindex1d(rz, z))
end

###

mutable struct Hist3_Cust
    edges::NTuple{3,Vector{CoordType}}
    nbins::NTuple{3,SizeType}
    origin::Vector3{CoordType}
    boxLength::Vector3{CoordType}
    weights::Array{SignalType,3}

    function Hist3_Cust(x::AbstractArray, y::AbstractArray, z::AbstractArray)
        nbins = (length(x) - 1, length(y) - 1, length(z) - 1)
        new(
            (x, y, z),
            nbins,
            V3[x[1], y[1], z[1]],
            V3[x[2] - x[1], y[2] - y[1], z[2] - z[1]],
            zeros(SignalType, nbins),
        )
    end
end

@inline function reset!(h::Hist3_Cust)
    h.weights = zeros(SignalType, h.nbins)
end

@propagate_inbounds function binindex1d(h::Hist3_Cust, d::SizeType, crd)
    dist = crd - h.origin[d]
    if dist < 0.0
        return 0
    end

    idx = floor(Int, dist / h.boxLength[d]) + 1
    if idx > h.nbins[d]
        return h.nbins[d] + 1
    end

    return idx
end

@propagate_inbounds function binindex(h::Hist3_Cust, x, y, z)
    return (binindex1d(h, 1, x), binindex1d(h, 2, y), binindex1d(h, 3, z))
end

@propagate_inbounds function atomic_push!(h::Hist3_Cust, x, y, z, wt)
    ix, iy, iz = binindex(h, x, y, z)
    lx, ly, lz = h.nbins
    if (unsigned(ix - 1) < lx) && (unsigned(iy - 1) < ly) && (unsigned(iz - 1) < lz)
        Atomix.@atomic h.weights[ix, iy, iz] += wt
    end
end

@inline binweights(h::Hist3_Cust) = h.weights

###

# const Hist3 = Hist3_FHist
const Hist3 = Hist3_Cust
