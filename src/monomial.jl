
#always assume the bases are sorted
struct monomial{N}
    coeff::Rational{Int64}
    factors::Array{Tuple{symmetric_polynomial{N},Int},1}
end

Base.:(==)(x::monomial{N},y::monomial{N}) where {N} = x.coeff == y.coeff && x.factors == y.factors

Base.convert(::Type{monomial{N}},x::symmetric_polynomial{N}) where {N} = monomial{N}(1,[(x,1)])
Base.convert(::Type{monomial{N}},x::elementary_symmetric_polynomial{N}) where {N} = monomial{N}(1,[(convert(symmetric_polynomial{N},x),1)])
function Base.convert(::Type{monomial{N}},x::elementary_monomial{N}) where {N}
    tmp = []
    for i =1:N
        if x.exponents[i]!=0
            push!(tmp,(convert(symmetric_polynomial{N},elementary_symmetric_polynomial{N}(i)),x.exponents[i]))
        end
    end
    return monomial{N}(x.coeff,tmp)
end

Base.promote_rule(::Type{monomial{N}},::Type{elementary_monomial{N}}) where {N} = monomial{N}
Base.promote_rule(::Type{monomial{N}},::Type{symmetric_polynomial{N}}) where {N} = monomial{N}
Base.promote_rule(::Type{monomial{N}},::Type{elementary_symmetric_polynomial{N}}) where {N} = monomial{N}
Base.promote_rule(::Type{elementary_monomial{N}},::Type{symmetric_polynomial{N}}) where {N} = monomial{N}

mergable(x::monomial{N},y::monomial{N}) where {N} = x.factors == y.factors
merge_unsafe(x::monomial{N},y::monomial{N}) where {N} = monomial{N}(x.coeff+y.coeff,x.factors)

function Base.:*(x::monomial{N},y::monomial{N}) where {N}
    x.factors==y.factors==[] && return monomial{N}(x.coeff*y.coeff,[])
    x.factors==[] && return monomial{N}(x.coeff*y.coeff,y.factors)
    y.factors==[] && return monomial{N}(x.coeff*y.coeff,x.factors)
    #use a double iterator method
    tmp_factors = []
    i,j=1,1
    while i <= length(x.factors) && j <= length(y.factors)
        if x.factors[i][1] < y.factors[j][1]
            push!(tmp_factors,x.factors[i])
            i+=1
        elseif y.factors[j][1] < x.factors[i][1]
            push!(tmp_factors,y.factors[j])
            j+=1
        else
            push!(tmp_factors,(x.factors[i][1],x.factors[i][2]+y.factors[j][2]))
            i+=1
            j+=1
        end
    end
    if i <= length(x.factors)
        append!(tmp_factors,x.factors[i:end])
    elseif j <= length(y.factors)
        append!(tmp_factors,y.factors[j:end])
    end
    return monomial{N}(x.coeff * y.coeff,tmp_factors)
end

Base.:*(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N} = convert(monomial{N},x) * convert(monomial{N},y)
Base.:*(x::symmetric_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::elementary_symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::symmetric_polynomial{N},y::elementary_monomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::elementary_monomial{N},y::symmetric_polynomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::monomial{N},y::symmetric_polynomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::monomial{N},y::elementary_monomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::monomial{N},y::elementary_symmetric_polynomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::symmetric_polynomial{N},y::monomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::elementary_monomial{N},y::monomial{N}) where {N} = *(promote(x,y)...)
Base.:*(x::elementary_symmetric_polynomial{N},y::monomial{N}) where {N} = *(promote(x,y)...)
Base.:*(c::T,x::monomial{N}) where {T<:Real} where {N} = monomial{N}(c*x.coeff,x.factors)
Base.:*(x::monomial{N},c::T) where {T<:Real} where {N} = monomial{N}(c*x.coeff,x.factors)
Base.:/(x::monomial{N},c::T) where {T<:Real} where {N} = monomial{N}(x.coeff//c,x.factors)
Base.:*(c::T,x::symmetric_polynomial{N}) where {T<:Real} where {N} = monomial{N}(c,[(x,1)])
Base.:*(x::symmetric_polynomial{N},c::T) where {T<:Real} where {N} = monomial{N}(c,[(x,1)])
Base.:/(x::symmetric_polynomial{N},c::T) where {T<:Real} where {N} = monomial{N}(1//c,[(x,1)])


Base.:-(x::symmetric_polynomial{N}) where {N} = monomial{N}(-1,[(x,1)])
Base.:-(x::monomial{N}) where {N} = monomial{N}(-x.coeff,x.factors)


function Base.isless(x::T,y::T) where T<:monomial
    for i = 1:length(x.factors)
        if x.factors[i] < y.factors[i]
            return true
        elseif x.factors[i] > y.factors[i]
            return false
        end
    end
    length(x.factors) < length(y.factors) && return true
    length(x.factors) > length(y.factors) && return false
    x.coeff < y.coeff && return true
    return false
end
