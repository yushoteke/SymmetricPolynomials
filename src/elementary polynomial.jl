
"""
An elementary polynomial is a polynomial with only elementary monomials,
which in turn is a monomial with only elementary symmetric polynomials.
"""
struct elementary_polynomial{N}
    terms::Array{elementary_monomial{N},1}
end
dim(x::elementary_polynomial{N}) where {N} = N

Base.convert(::Type{elementary_polynomial{N}},x::elementary_monomial{N}) where {N} = elementary_polynomial{N}([x])
Base.convert(::Type{elementary_polynomial{N}},x::elementary_symmetric_polynomial{N}) where {N} = elementary_polynomial{N}([convert(elementary_monomial{N},x)])
Base.promote_rule(::Type{elementary_polynomial{N}},::Type{elementary_monomial{N}}) where {N} = elementary_polynomial{N}
Base.promote_rule(::Type{elementary_polynomial{N}},::Type{elementary_symmetric_polynomial{N}}) where {N} = elementary_polynomial{N}

Base.:(==)(x::elementary_polynomial{N},y::elementary_polynomial{N}) where {N} = x.terms == y.terms

function simplify_monomial_array(arr)
    if length(arr) <= 1
        arr[end].coeff == 0 && return eltype(arr)[]
        return arr
    end
    tmp = sort(arr)
    result = eltype(arr)[tmp[1]]
    for i in tmp[2:end]
        if mergable(i,result[end])
            result[end] = merge_unsafe(result[end],i)
        elseif result[end].coeff == 0
            result[end] = i
        else
            push!(result,i)
        end
    end
    result[end].coeff == 0 && pop!(result)
    return result
end

#ep + ep = ep
function Base.:+(x::elementary_polynomial{N},y::elementary_polynomial{N}) where {N}
    length(x.terms)==0 && return y
    length(y.terms)==0 && return x

    tmp = eltype(x.terms)[]
    i,j=1,1
    while i<=length(x.terms) && j<=length(y.terms)
        if mergable(x.terms[i],y.terms[j])
            x.terms[i].coeff != -y.terms[j].coeff && push!(tmp,merge_unsafe(x.terms[i],y.terms[j]))
            i,j=i+1,j+1
        elseif x.terms[i] < y.terms[j]
            push!(tmp,x.terms[i])
            i+=1
        elseif x.terms[i] > y.terms[j]
            push!(tmp,y.terms[j])
            j+=1
        end
    end
    i <=length(x.terms) && append!(tmp,x.terms[i:end])
    j <=length(y.terms) && append!(tmp,y.terms[j:end])
    return elementary_polynomial{N}(tmp)
end
Base.:+(x::elementary_symmetric_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = convert(elementary_polynomial{N},x) + convert(elementary_polynomial{N},y)
Base.:+(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = convert(elementary_polynomial{N},x) + convert(elementary_polynomial{N},y)

#handles esp + ep = ep and em + ep = ep
Base.:+(x::elementary_polynomial{N},y::elementary_monomial{N}) where {N} = x + convert(elementary_polynomial{N},y)
Base.:+(y::elementary_monomial{N},x::elementary_polynomial{N}) where {N} = x + convert(elementary_polynomial{N},y)
Base.:+(x::elementary_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = x + convert(elementary_polynomial{N},y)
Base.:+(y::elementary_symmetric_polynomial{N},x::elementary_polynomial{N}) where {N} = x + convert(elementary_polynomial{N},y)
Base.:+(x::elementary_monomial{N},y::elementary_symmetric_polynomial{N}) where {N} = x + convert(elementary_monomial{N},y)
Base.:+(y::elementary_symmetric_polynomial{N},x::elementary_monomial{N}) where {N} = x + convert(elementary_polynomial{N},y)

function Base.:-(x::elementary_polynomial{N},y::elementary_polynomial{N}) where {N}
    length(x.terms)==0 && return y
    length(y.terms)==0 && return x

    tmp = eltype(x.terms)[]
    i,j=1,1
    while i<=length(x.terms) && j<=length(y.terms)
        if mergable(x.terms[i],y.terms[j])
            x.terms[i].coeff != y.terms[j].coeff && push!(tmp,merge_unsafe(x.terms[i],-y.terms[j]))
            i,j=i+1,j+1
        elseif x.terms[i] < y.terms[j]
            push!(tmp,x.terms[i])
            i+=1
        elseif x.terms[i] > y.terms[j]
            push!(tmp,-y.terms[j])
            j+=1
        end
    end
    i <=length(x.terms) && append!(tmp,x.terms[i:end])
    j <=length(y.terms) && append!(tmp,.-y.terms[j:end])
    return elementary_polynomial{N}(tmp)
end
Base.:-(x::elementary_polynomial{N}) where {N} = elementary_polynomial{N}(.-x.terms)
Base.:-(x::elementary_symmetric_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = convert(elementary_polynomial{N},x) + (-1*y)
Base.:-(x::elementary_monomial{N},y::elementary_monomial{N}) where {N} = convert(elementary_polynomial{N},x) + (-1*y)

Base.:-(x::elementary_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = x - convert(elementary_polynomial{N},y)
Base.:-(x::elementary_symmetric_polynomial{N},y::elementary_polynomial{N}) where {N} = convert(elementary_polynomial{N},x) - y
Base.:-(x::elementary_polynomial{N},y::elementary_monomial{N}) where {N} = x - convert(elementary_polynomial{N},y)
Base.:-(x::elementary_monomial{N},y::elementary_polynomial{N}) where {N} = convert(elementary_polynomial{N},x) - y
Base.:-(x::elementary_monomial{N},y::elementary_symmetric_polynomial{N}) where {N} = convert(elementary_polynomial{N},x) - y
Base.:-(x::elementary_symmetric_polynomial{N},y::elementary_monomial{N}) where {N} = convert(elementary_polynomial{N},x) - y

"""
    multiply_karatsubalike(x,y)

For the naive implementation, say x.terms has length n and y.terms has length m,
Then the total runtime is n*m*log(nm),with the quicksort at the end

However, we can also take advantage of the identity, say x=x1+x2,y=y1+y2
Then u=(x1+y1)(x2+y2)=x1x2+x1y2+y1x2+y1y2,v=(x1-y2)(x2-y1)=x1x2-x1y1-x2y2+y1y2
u-v=x1y2+y1x2+x1y1+x2y2=xy
Lets say we divide x into two equal halves, so that length of x1 is n/2, and similarly for
y, length of y1 is m/2, the temporary result x1+y1 has at most length n/2+m/2, but due to simplifications,
It's highly likely that the length is smaller than n/2+m/2.
For example, when we decomposed X=symmetric_polynomial(10,9...1,0), the polynomial has 2485 terms, and
decomposed Y=symmetric_polynomial(11,8,8,7,6...1,0) has 3023 terms. But if we sum up the first 1300 terms
of X and the first 1500 terms of Y, the resulting polynomial has only 1500 terms! This means they completely
merged together. What's even more amazing, summing the terms 1301:end of X and 1501 of Y, the resulting
polynomial even got shorter. This means some terms completely cancelled out.



"""

Base.:*(y::elementary_polynomial{N},x::elementary_monomial{N}) where {N} = elementary_polynomial{N}(simplify_monomial_array([x*i for i in y.terms]))
Base.:*(x::elementary_monomial{N},y::elementary_polynomial{N}) where {N} = y * x
Base.:*(x::elementary_polynomial{N},y::elementary_polynomial{N}) where {N} = sum([i*y for i in x.terms])

function Base.:^(x::elementary_polynomial{N},y::Int64) where {N}
    y<1 && throw(DomainError("power argument should >=1"))
    if y==1
        return x
    elseif y%2 == 0
        return (x*x)^(y÷2)
    else
        return x * (x^(y-1))
    end
end

Base.:*(c::T,y::elementary_polynomial{N}) where {T<:Real,N} = elementary_polynomial{N}(c.*y.terms)
Base.:*(y::elementary_polynomial{N},c::T) where {T<:Real,N} = elementary_polynomial{N}(c.*y.terms)
Base.:/(y::elementary_polynomial{N},c::T) where {T<:Real,N} = elementary_polynomial{N}((1//c).*y.terms)
Base.:*(x::elementary_polynomial{N},y::elementary_symmetric_polynomial{N}) where {N} = elementary_polynomial{N}(simplify_monomial_array([y*i for i in x.terms]))
Base.:*(y::elementary_symmetric_polynomial{N},x::elementary_polynomial{N}) where {N} = elementary_polynomial{N}(simplify_monomial_array([y*i for i in x.terms]))
