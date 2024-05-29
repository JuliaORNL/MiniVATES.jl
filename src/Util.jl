import StaticArrays
import StaticArrays: SVector, MVector
import JACC
import Base: @propagate_inbounds

@inline function swap(a, b)
    return (b, a)
end

@inline function lerp(a, b, t)
    return a + t * (b - a)
end

const SignalType = Float32
const CoordType = Float32
const ScalarType = Float32

const Vector3{T} = StaticArrays.SVector{3,T} where {T}
const Vector4{T} = StaticArrays.SVector{4,T} where {T}

const Vec3 = Vector3{ScalarType}
@inline Vec3() = Vec3(0, 0, 0)
const V3 = StaticArrays.SA{ScalarType}
const Vec4 = Vector4{ScalarType}
const V4 = StaticArrays.SA{ScalarType}

const Crd3 = Vector3{CoordType}
@inline Crd3() = Crd3(0, 0, 0)
const Crd4 = Vector4{CoordType}
const C3 = StaticArrays.SA{CoordType}
const C4 = StaticArrays.SA{CoordType}

# const SizeType = UInt64
const SizeType = Int
const Id3 = StaticArrays.SVector{3,SizeType}
@inline Id3() = Vec3(0, 0, 0)
const I3 = StaticArrays.SA{SizeType}

function setUpIndexMaker(indexMax::Id3)
    out = ones(StaticArrays.MVector{3,SizeType})
    for i = 2:3
        out[i] = out[i - 1] * indexMax[i - 1]
    end
    return out
end

const SquareMatrix3{T} = StaticArrays.SMatrix{3,3,T} where {T}
@inline SquareMatrix3{T}() where {T} = zeros(SquareMatrix3{T})

const SquareMatrix3V = SquareMatrix3{ScalarType}
const SquareMatrix3C = SquareMatrix3{CoordType}

const Array1 = JACC.Array{T,1} where {T}
@inline Array1{T}(n::Int64) where {T} = Array1{T}(undef, n)
const Array1r = Array1{ScalarType}

mutable struct PreallocVector{T} <: AbstractVector{T}
    data::Vector{T}
    # maxIdx::SizeType
    curIdx::SizeType

    PreallocVector(v::Vector{T}) where {T} = new{T}(v, 1)
end

@inline data(v::PreallocVector) = v.data

@inline Base.convert(::Type{PreallocVector}, v) = PreallocVector(v)

@inline Base.length(v::PreallocVector) = v.curIdx - 1
@inline Base.size(v::PreallocVector) = (length(v))

@propagate_inbounds Base.getindex(v::PreallocVector, i) = v.data[i]

@propagate_inbounds function Base.setindex!(v::PreallocVector, x, i)
    setindex!(v.data, x, i)
end

@propagate_inbounds function Base.push!(v::PreallocVector, val)
    v.data[v.curIdx] = val
    v.curIdx += 1
    return v
end

@inline function Base.fill!(v::PreallocVector, val, len)
    fill!(@view(v.data[1:len]), val)
    return resize!(v, len)
end

@inline function Base.resize!(v::PreallocVector, len)
    v.curIdx = len + 1
    return v
end

@inline function Base.empty!(v::PreallocVector)
    v.curIdx = 1
end

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
