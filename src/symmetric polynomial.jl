
struct symmetric_polynomial{N}
    exponents::NTuple{N,Int64}
end

"""
    symmetric_polynomial(x...)

Examples:
symmetric_polynomial(1,0,0)     = x + y + z
symmetric_polynomial(1,0,0,0)   = w + x + y + z
symmetric_polynomial(1,1,0)     = xy + xz + yz
symmetric_polynomial(1,1,1)     = xyz
symmetric_polynomial(2,1,0)     = x^2y + x^2z + y^2x + y^2z + z^2x + z^2y
"""

function symmetric_polynomial(x...)
    tmp = length(x)<10 ? TupleTools.sort(x) : tuple(sort([x...])...)
    return symmetric_polynomial{length(x)}(tmp)
end

is_elementary(x::T) where T<:symmetric_polynomial = x.exponents[end]<=1
dim(x::symmetric_polynomial{N}) where {N} = N
Base.:<(x::T,y::T) where T<:symmetric_polynomial = x.exponents < y.exponents
Base.:>(x::T,y::T) where T<:symmetric_polynomial = x.exponents > y.exponents
Base.:(<=)(x::T,y::T) where T<:symmetric_polynomial = x.exponents <= y.exponents
Base.:(>=)(x::T,y::T) where T<:symmetric_polynomial = x.exponents >= y.exponents
Base.isless(x::T,y::T) where T<:symmetric_polynomial = x.exponents < y.exponents
