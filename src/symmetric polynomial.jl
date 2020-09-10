
struct symmetric_polynomial{N}
    exponents::NTuple{N,Int64}
end

function symmetric_polynomial(x...)
    tmp = length(x)<10 ? TupleTools.sort(x) : tuple(sort([x...])...)
    return symmetric_polynomial{length(x)}(tmp)
end

Base.convert(::Type{symmetric_polynomial{N}},x::elementary_symmetric_polynomial{N}) where {N} = symmetric_polynomial{N}(ntuple(i->i<=N-x.order ? 0 : 1,N))
Base.promote_rule(::Type{elementary_symmetric_polynomial{N}},::Type{symmetric_polynomial{N}}) where {N} = symmetric_polynomial{N}
is_elementary(x::symmetric_polynomial{N}) where {N} = x.exponents[end]<=1
dim(x::symmetric_polynomial{N}) where {N} = N
Base.isless(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N} = x.exponents < y.exponents
mergable(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N} = x.exponents == y.exponents
"""
    num_highest(x::symmetric_polynomial{N})

Find the number of highest elements.
For example, S(0,1,1,2) = 1, S(0,1,1,1) = 3, S(0,1,2,2) = 2

"""
num_highest(x::symmetric_polynomial{N}) where {N} = length(findall(i->i==x.exponents[end],x.exponents))
