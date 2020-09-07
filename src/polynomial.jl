
struct polynomial{N}
    terms::Array{monomial{N},1}
end

polynomial(x::Array{T,1}) where {T<:monomial} = polynomial{N}(x)
to_polynomial(x::T) where {T<:monomial} = polynomial([x])
to_polynomial(x::T) where {T<:symmetric_polynomial} = to_polynomial(to_monomial(x))

function simplify_monomial_array(arr::Array{T,1}) where T<:monomial
    length(arr) <= 1 && return arr
    tmp = sort(arr)
    result = [tmp[1]]
    for i in tmp[2:end]
        if summable(i,result[end])
            result[end] += i
        else
            push!(result,i)
        end
    end
    return result[findall(t->t.coeff!=0,result)]
end

Base.:+(x::T,y::T) where {T<:polynomial} = polynomial(simplify_monomial_array(cat(x.terms,y.terms,dims=1)))
Base.:-(x::T) where {T<:polynomial} = polynomial(.-x.terms)
Base.:-(x::T,y::T) where {T<:polynomial} = polynomial(simplify_monomial_array(cat(x.terms,.-y.terms,dims=1)))

function Base.:*(x::T,y::T) where {T<:polynomial}
    tmp = eltype(x.terms)[]
    for i in x.terms
        for j in y.terms
            push!(tmp,i*j)
        end
    end
    return polynomial(simplify_monomial_array(tmp))
end

function Base.:^(x::polynomial,y::Int64)
    y<1 && throw(DomainError("power argument should >=1"))
    if y==1
        return x
    elseif y%2 == 0
        return (x*x)^(y÷2)
    else
        return x * (x^(y-1))
    end
end

Base.:*(y::polynomial{N},x::monomial{N}) where {N} = y * to_polynomial(x)
Base.:*(x::monomial{N},y::polynomial{N}) where {N} = y * to_polynomial(x)

Base.:*(c::T,y::S) where {T<:Real,S<:polynomial} = polynomial(c.*y.terms)
Base.:*(y::S,c::T) where {T<:Real,S<:polynomial} = polynomial(c.*y.terms)
