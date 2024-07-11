import Adapt
import Base: @propagate_inbounds

struct PreallocJaggedArray{T,TDataVec,TIdVec}
    data::TDataVec
    rowSize::SizeType
    start::TIdVec
    curIdx::TIdVec
end

function PreallocJaggedArray{T}(o::TV1, rowSize::SizeType, s::TV2, i::TV2) where {T,TV1,TV2}
    PreallocJaggedArray{T,TV1,TV2}(o, rowSize, s, i)
end

function Adapt.adapt_structure(to, a::PreallocJaggedArray{T}) where {T}
    PreallocJaggedArray{T}(
        adapt_structure(to, a.data),
        a.rowSize,
        adapt_structure(to, a.start),
        adapt_structure(to, a.curIdx),
    )
end

function PreallocJaggedArray{T}(
    m::Array1{T},
    rowCount::SizeType,
    rowSize::SizeType,
) where {T}
    PreallocJaggedArray{T,Array1{T},Array1{SizeType}}(
        m,
        rowSize,
        Array1([((i - 1) * rowSize) for i = 1:rowCount]),
        JACC.ones(SizeType, rowCount),
    )
end

@inline function PreallocJaggedArray{T}(rowCount, rowSize) where {T}
    PreallocJaggedArray{T}(Array1{T}(undef, rowCount * rowSize), rowCount, rowSize)
end

@inline function PreallocJaggedArray{T}() where {T}
    PreallocJaggedArray{T}(Array1{T}(), 0, 0, 0)
end

function reset!(a::PreallocJaggedArray)
    if a.rowSize == 0
        return nothing
    end
    unsafe_free!(a.data)
    unsafe_free!(a.start)
    unsafe_free!(a.curIdx)
    return nothing
end

@propagate_inbounds function rowview(a::PreallocJaggedArray, n)
    start = a.start[n] + 1
    stop = a.start[n] + length(a, n)
    return view(a.data, start:stop)
end

@inline Base.length(a::PreallocJaggedArray) = length(a.start)

@inline rowSize(a::PreallocJaggedArray) = a.rowSize

@propagate_inbounds function Base.length(a::PreallocJaggedArray, n)
    a.curIdx[n] - 1
end
@propagate_inbounds function Base.getindex(a::PreallocJaggedArray, n, i)
    getindex(a.data, a.start[n] + i)
end
@propagate_inbounds function Base.setindex!(a::PreallocJaggedArray, x, n, i)
    setindex!(a.data, x, a.start[n] + i)
end

@propagate_inbounds function Base.push!(a::PreallocJaggedArray, n, val)
    a.data[a.start[n] + a.curIdx[n]] = val
    a.curIdx[n] += 1
    return a
end

@propagate_inbounds function Base.fill!(a::PreallocJaggedArray, n, val, len)
    resize!(a, n, len)
    fill!(rowview(a, n), val)
    return a
end

@propagate_inbounds function Base.resize!(a::PreallocJaggedArray, n, len)
    a.curIdx[n] = len + 1
    return a
end

@propagate_inbounds function Base.empty!(a::PreallocJaggedArray, n)
    a.curIdx[n] = 1
    return a
end

struct PreallocArrayRow{T,TV1,TV2} <: AbstractVector{T}
    a::PreallocJaggedArray{T,TV1,TV2}
    n::SizeType
end

@propagate_inbounds function row(a::PreallocJaggedArray{T,TV1,TV2}, n) where {T,TV1,TV2}
    return PreallocArrayRow{T,TV1,TV2}(a, n)
end

@propagate_inbounds Base.length(r::PreallocArrayRow) = length(r.a, r.n)
@propagate_inbounds Base.size(r::PreallocArrayRow) = (length(r),)
@propagate_inbounds Base.firstindex(r::PreallocArrayRow) = 1
@propagate_inbounds Base.lastindex(r::PreallocArrayRow) = length(r)
@propagate_inbounds Base.getindex(r::PreallocArrayRow, i) = getindex(r.a, r.n, i)
@propagate_inbounds function Base.setindex!(r::PreallocArrayRow, x, i)
    setindex!(r.a, x, r.n, i)
    return r
end
@propagate_inbounds function Base.push!(r::PreallocArrayRow, val)
    push!(r.a, r.n, val)
    return r
end
@propagate_inbounds function Base.empty!(r::PreallocArrayRow)
    empty!(r.a, r.n)
    return r
end
@propagate_inbounds function Base.fill!(r::PreallocArrayRow, val, len)
    fill!(r.a, r.n, val, len)
    return r
end
@propagate_inbounds function Base.resize!(r::PreallocArrayRow, len)
    resize!(r.a, r.n, len)
    return r
end

# @inline function Base.sortperm!(ix::PreallocArrayRow, v::PreallocArrayRow; lt)
#     len = length(v)
#     if len > 0
#         if length(ix) != len
#             resize!(ix, len)
#         end
#         @inbounds sortperm!(
#             rowview(ix.a, ix.n),
#             rowview(v.a, v.n),
#             lt = lt,
#             alg = MergeSort,
#         )
#     end
#     return (ix, v)
# end

struct SortedPreallocRow{T} <: AbstractVector{T}
    perm::PreallocArrayRow{SizeType}
    data::PreallocArrayRow{T}
end


# struct SortedPreallocVector{T} <: AbstractVector{T}
#     perm::PreallocVector{SizeType}
#     data::PreallocVector{T}

# function SortedPreallocVector(
#     p::PreallocVector{SizeType},
#     d::PreallocVector{T},
# ) where {T}
#     new{T}(sortperm!(p, d))
# end
# end

# @propagate_inbounds function Base.getindex(v::SortedPreallocVector, i::Integer)
#     return v.data[v.perm[i]]
# end

# @inline Base.length(v::SortedPreallocVector) = length(v.data)

