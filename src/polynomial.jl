
struct polynomial{N}
    terms::Array{monomial{N},1}
end

dim(x::polynomial{N}) where {N} = N
polynomial(x::Array{monomial{N},1}) where {N} = polynomial{N}(x)
to_polynomial(x::monomial{N}) where {N} = polynomial(monomial{N}[x])
to_polynomial(x::symmetric_polynomial{N}) where {N} = to_polynomial(to_monomial(x))

function simplify_monomial_array(arr::Array{monomial{N},1}) where {N}
    length(arr) <= 1 && return arr
    tmp = sort(arr)
    result = [tmp[1]]
    for i in tmp[2:end]
        if summable(i,result[end])
            result[end] += i
        elseif result[end].coeff == 0
            result[end] = i
        else
            push!(result,i)
        end
    end
    result[end].coeff == 0 && pop!(result)
    return result
end

function Base.:+(x::polynomial{N},y::polynomial{N}) where {N}
    tmp = copy(x.terms)
    append!(tmp,y.terms)
    return polynomial(simplify_monomial_array(tmp))
end
Base.:-(x::polynomial{N}) where {N} = polynomial(.-x.terms)
Base.:-(x::polynomial{N},y::polynomial{N}) where {N} = x + (-y)

Base.:*(y::polynomial{N},x::monomial{N}) where {N} = polynomial(simplify_monomial_array([x*i for i in y.terms]))
Base.:*(x::monomial{N},y::polynomial{N}) where {N} = y * x
Base.:*(x::polynomial{N},y::polynomial{N}) where {N} = sum([i*y for i in x.terms])

function Base.:^(x::polynomial{N},y::Int64) where {N}
    y<1 && throw(DomainError("power argument should >=1"))
    if y==1
        return x
    elseif y%2 == 0
        return (x*x)^(y÷2)
    else
        return x * (x^(y-1))
    end
end

Base.:*(c::T,y::polynomial{N}) where {T<:Real,N} = polynomial(c.*y.terms)
Base.:*(y::polynomial{N},c::T) where {T<:Real,N} = polynomial(c.*y.terms)
