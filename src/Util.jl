import StaticArrays
import StaticArrays: SVector, MVector
import JACC

@inline function swap(a, b)
    return (b, a)
end

@inline function lerp(a, b, t)
    return a + t * (b - a)
end

# const SignalType = Float64

const ScalarType = Float64
const Vec3 = StaticArrays.SVector{3,ScalarType}
@inline Vec3() = Vec3(0, 0, 0)
const V3 = StaticArrays.SA{ScalarType}
const Vec4 = StaticArrays.SVector{4,ScalarType}
const V4 = StaticArrays.SA{ScalarType}

const CoordType = Float32
const Crd3 = StaticArrays.SVector{3,CoordType}
@inline Crd3() = Crd3(0, 0, 0)
const C3 = StaticArrays.SA{CoordType}

# const SizeType = UInt64
const SizeType = Int
const Id3 = StaticArrays.SVector{3,SizeType}
@inline Id3() = Vec3(0, 0, 0)
const I3 = StaticArrays.SA{SizeType}

function setUpIndexMaker(indexMax::Id3)
    out = ones(StaticArrays.MVector{3,SizeType})
    for i = 2:3
        out[i] = out[i-1] * indexMax[i-1]
    end
    return out
end

const SquareMatrix3 = StaticArrays.SMatrix{3,3,ScalarType}
# const SquareMatrix3 = StaticArrays.SMatrix{Tuple{3,3},ScalarType,2,9}
@inline SquareMatrix3() = zeros(SquareMatrix3)

const Array1 = JACC.Array{T,1} where {T}
@inline Array1{T}(n::Int64) where {T} = Array1{T}(undef, n)
const Array1r = Array1{ScalarType}

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
