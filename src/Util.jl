import StaticArrays
import StaticArrays: SVector, MVector
import JACC
import Adapt
import Base: @propagate_inbounds

@inline function lerp(a, b, t)
    return a + t * (b - a)
end

const SignalType = Float32
const CoordType = Float32
const ScalarType = Float32
const SizeType = Int

const Vector3{T} = StaticArrays.SVector{3,T} where {T}
@inline Vector3{T}() where {T} = Vector3{T}(0, 0, 0)
const Vector4{T} = StaticArrays.SVector{4,T} where {T}

const Vec3 = Vector3{ScalarType}
const V3 = StaticArrays.SA{ScalarType}
const Vec4 = Vector4{ScalarType}
const V4 = StaticArrays.SA{ScalarType}

const Crd3 = Vector3{CoordType}
const Crd4 = Vector4{CoordType}
const C3 = StaticArrays.SA{CoordType}
const C4 = StaticArrays.SA{CoordType}

const Id3 = Vector3{SizeType}
const I3 = StaticArrays.SA{SizeType}

function setUpIndexMaker(indexMax::Id3)
    out = ones(StaticArrays.MVector{3,SizeType})
    for i = 2:3
        out[i] = out[i - 1] * indexMax[i - 1]
    end
    return out
end

# const SquareMatrix3{T} = StaticArrays.SMatrix{3,3,T} where {T}
const SquareMatrix3{T} = StaticArrays.SMatrix{3,3,T,9} where {T}
@inline SquareMatrix3{T}() where {T} = zeros(SquareMatrix3{T})

const SquareMatrix3r = SquareMatrix3{ScalarType}
const SquareMatrix3c = SquareMatrix3{CoordType}

const Array1 = JACC.Array{T,1} where {T}
@inline Array1{T}(n::Int64) where {T} = Array1{T}(undef, n)
const Array1r = Array1{ScalarType}
const Array1c = Array1{CoordType}

const Array2 = JACC.Array{T,2} where {T}
@inline Array2{T}(m::Int64, n::Int64) where {T} = Array2{T}(undef, m, n)
const Array2c = Array2{CoordType}

using Core: LLVMPtr

struct PreallocJaggedArray{T, TDataVec, TIdVec}
    data::TDataVec
    rowSize::SizeType
    start::TIdVec
    curIdx::TIdVec
end

function PreallocJaggedArray{T}(o::TV1, rowSize::SizeType, s::TV2, i::TV2) where {T, TV1, TV2}
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

function PreallocJaggedArray{T}(m::Array1{T}, rowCount::SizeType, rowSize::SizeType) where {T}
    PreallocJaggedArray{T, Array1{T}, Array1{SizeType}}(
        m,
        rowSize,
        Array1([((i - 1) * rowSize) for i = 1:rowCount]),
        JACC.ones(SizeType, rowCount),
    )
end

@inline function PreallocJaggedArray{T}(rowCount, rowSize) where {T}
    PreallocJaggedArray{T}(Array1{T}(undef, rowCount * rowSize), rowCount, rowSize)
end

@propagate_inbounds function rowview(a::PreallocJaggedArray, n)
    start = a.start[n] + 1
    stop = length(a, n)
    return view(a.data, start:stop)
end

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
    a.data[a.curIdx[n]] = val
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

@inline function Base.sortperm!(ix::PreallocArrayRow, v::PreallocArrayRow; lt)
    len = length(v)
    if len > 0
        if length(ix) != len
            resize!(ix, len)
        end
        @inbounds sortperm!(rowview(ix.a, ix.n), rowview(v.a, v.n), lt = lt,
                           alg = MergeSort)
    end
    return (ix, v)
end

struct SortedPreallocRow{T} <: AbstractVector{T}
    perm::PreallocArrayRow{SizeType}
    data::PreallocArrayRow{T}
end

# mutable struct PreallocVector{T} <: CuDeviceVector{T}
struct PreallocVector{T} <: DenseMatrix{T}
    data::Array1{T}
    curIdx::Array1{SizeType}

    PreallocVector(v::Array1{T}) where {T} = new{T}(v, 1)

    function PreallocVector(v::CUDA.CuDeviceVector{T}, row, rowSize) where {T}
        start = (row - 1) * rowSize + 1
        # stop = start + rowSize - 1
        new{T}(pointer(v, start), 1)
    end
end

@inline data(v::PreallocVector) = v.data

@inline Base.convert(::Type{PreallocVector}, v) = PreallocVector(v)

@inline Base.length(v::PreallocVector) = v.curIdx - 1
@inline Base.size(v::PreallocVector) = (length(v),)

@propagate_inbounds Base.getindex(v::PreallocVector, i) = unsafe_load(v.data, i)

@propagate_inbounds function Base.setindex!(v::PreallocVector, x, i)
    unsafe_store!(v.data, x, i)
end

# @propagate_inbounds function Base.push!(v::PreallocVector, val)
#     v.data[v.curIdx] = val
#     v.curIdx += 1
#     return v
# end

# @inline function Base.fill!(v::PreallocVector, val, len)
#     fill!(@view(v.data[1:len]), val)
#     return resize!(v, len)
# end

# @inline function Base.resize!(v::PreallocVector, len)
#     v.curIdx = len + 1
#     return v
# end

# @inline function Base.empty!(v::PreallocVector)
#     v.curIdx = 1
# end

@inline function Base.sortperm!(ix::PreallocVector, v::PreallocVector; lt)
    len = length(v)
    if len > 0
        @inbounds sortperm!(@view(ix.data[1:len]), @view(v.data[1:len]), lt = lt)
    end
    return (ix, v)
end

struct SortedPreallocVector{T} <: AbstractVector{T}
    perm::PreallocVector{SizeType}
    data::PreallocVector{T}

    # function SortedPreallocVector(
    #     p::PreallocVector{SizeType},
    #     d::PreallocVector{T},
    # ) where {T}
    #     new{T}(sortperm!(p, d))
    # end
end

@propagate_inbounds function Base.getindex(v::SortedPreallocVector, i::Integer)
    return v.data[v.perm[i]]
end

@inline Base.length(v::SortedPreallocVector) = length(v.data)

struct Optional{T}
    value::Union{T,Missing}
    Optional{T}() where {T} = new(missing)
    Optional{T}(x) where {T} = new(x)
    Optional{T}(::Missing) where {T} = new(missing)
end

ismissing(x::Optional) = x.value === missing
hasvalue(x::Optional) = !ismissing(x)
coalesce(x::Optional) = x.value
Base.convert(::Type{Optional{T}}, ::Missing) where {T} = Optional{T}()
