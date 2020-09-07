
struct symmetric_polynomial{N}
    exponents::NTuple{N,Int64}
end

function symmetric_polynomial(x...)
    tmp = length(x)<10 ? TupleTools.sort(x) : tuple(sort([x...])...)
    return symmetric_polynomial{length(x)}(tmp)
end

symmetric_polynomial(x::Array{Int,1}) = symmetric_polynomial{length(x)}(tuple(sort(x)...))

is_elementary(x::T) where T<:symmetric_polynomial = x.exponents[end]<=1
dim(x::symmetric_polynomial{N}) where {N} = N
Base.:<(x::T,y::T) where T<:symmetric_polynomial = x.exponents < y.exponents
Base.:>(x::T,y::T) where T<:symmetric_polynomial = x.exponents > y.exponents
Base.:(<=)(x::T,y::T) where T<:symmetric_polynomial = x.exponents <= y.exponents
Base.:(>=)(x::T,y::T) where T<:symmetric_polynomial = x.exponents >= y.exponents
Base.isless(x::T,y::T) where T<:symmetric_polynomial = x.exponents < y.exponents
