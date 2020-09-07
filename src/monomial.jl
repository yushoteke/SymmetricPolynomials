
#always assume the bases are sorted
struct monomial{N}
    coeff::Rational{Int64}
    factors::Tuple{Vararg{Tuple{symmetric_polynomial{N},Int}}}
end

monomial(c,x::symmetric_polynomial{N}) where {N} = monomial{N}(c,((x,1),))
summable(x::T,y::T) where T<:monomial = x.factors == y.factors
to_monomial(x::symmetric_polynomial{N}) where {N} = monomial(1//1,x)
Base.:+(x::monomial{N},y::monomial{N}) where {N} =
    summable(x,y) ? monomial{N}(x.coeff + y.coeff,x.factors) : throw(DomainError("input factors need to be same"))

function Base.:*(x::monomial{N},y::monomial{N}) where {N}
    x.factors==y.factors==tuple() && return monomial{N}(x.coeff*y.coeff,tuple())
    x.factors==tuple() && return monomial{N}(x.coeff*y.coeff,y.factors)
    y.factors==tuple() && return monomial{N}(x.coeff*y.coeff,x.factors)
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
        push!(tmp_factors,x.factors[i:end]...)
    elseif j <= length(y.factors)
        push!(tmp_factors,y.factors[j:end]...)
    end
    return monomial{N}(x.coeff * y.coeff,tuple(tmp_factors...))
end

Base.:*(x::T,y) where T<:monomial = monomial(x.coeff * y,x.factors)
Base.:*(y,x::T) where T<:monomial = monomial(x.coeff * y,x.factors)
Base.:-(x::T) where T<:monomial = monomial(-x.coeff,x.factors)

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
