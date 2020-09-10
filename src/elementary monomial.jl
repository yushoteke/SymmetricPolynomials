
"""
elementary monomials coeff,(a,b,c...) represents
coeff * e1^a * e2^b * e3^c ...
"""
struct elementary_monomial{N}
    coeff::Rational{Int64}
    exponents::NTuple{N,Int64}
end
dim(x::elementary_monomial{N}) where {N} = N

Base.convert(::Type{elementary_monomial{N}},x::elementary_symmetric_polynomial{N}) where {N} = elementary_monomial{N}(1,ntuple(i->i==x.order ? 1 : 0,N))
Base.promote_rule(::Type{elementary_monomial{N}},::Type{elementary_symmetric_polynomial{N}}) where {N} = elementary_monomial{N}

mergable(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = x.exponents == y.exponents
merge_unsafe(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = elementary_monomial{N}(x.coeff + y.coeff,x.exponents)

#esp * esp = em
Base.:*(x::elementary_symmetric_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = convert(elementary_monomial{N},x) * convert(elementary_monomial{N},y)
Base.:*(x::elementary_symmetric_polynomial{N},y::elementary_monomial{N}) where {N} = y * convert(elementary_monomial{N},x)
Base.:*(y::elementary_monomial{N},x::elementary_symmetric_polynomial{N}) where {N} = y * convert(elementary_monomial{N},x)
#em * em = em
Base.:*(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = elementary_monomial{N}(x.coeff * y.coeff,x.exponents .+ y.exponents)
#em * c = em
Base.:*(x::elementary_monomial{N},y::T) where {N,T<:Real} = elementary_monomial{N}(x.coeff * y,x.exponents)
Base.:*(y::T,x::elementary_monomial{N}) where {N,T<:Real} = x * y
Base.:/(x::elementary_monomial{N},y::T) where {N,T<:Real} = x * (1//y)
Base.:-(x::elementary_monomial{N}) where {N} = x * (-1)
#esp * c = em
Base.:*(x::elementary_symmetric_polynomial{N},y::T) where {N,T<:Real} = convert(elementary_monomial{N},x) * y
Base.:*(y::T,x::elementary_symmetric_polynomial{N}) where {N,T<:Real} = convert(elementary_monomial{N},x) * y
Base.:/(x::elementary_symmetric_polynomial{N},y::T) where {N,T<:Real} = x * (1//y)
Base.:-(x::elementary_symmetric_polynomial{N}) where {N} = x * (-1)


Base.isless(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = x.exponents != y.exponents ? x.exponents < y.exponents : x.coeff < y.coeff
