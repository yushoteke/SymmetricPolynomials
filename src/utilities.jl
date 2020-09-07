
function ordering_m1_n0(m,n)
    #generates all unique orderings or m 1s and n 0s
    #=
    For Example, 2 1s and 2 0s
    1100,1010,1001,0110,0101,0011
    =#
    #base case
    if m<=n<=0
        return []
    elseif m<=0 && n>=1
        return [zeros(Int64,n)]
    elseif m>=1 && n<=0
        return [ones(Int64,m)]
    else    #both m,n>=1
        #assume first entry is 0
        tmp0 = ordering_m1_n0(m,n-1)
        for i in tmp0
            pushfirst!(i,0)
        end
        #assume first entry is 1
        tmp1 = ordering_m1_n0(m-1,n)
        for i in tmp1
            pushfirst!(i,1)
        end
        append!(tmp0,tmp1)
        return tmp0
    end
end

function get_valid_exponents(dim,rank)
    tmp = integer_partitions(rank)
    too_long = findall(x->length(x)>dim,tmp)
    deleteat!(tmp,too_long)
    for i in tmp
        if length(i) < dim
            append!(i,zeros(Int64,dim-length(i)))
        end
        reverse!(i)
    end
    return tmp
end

function symmetry_factor(exponents)
    tmp = counter(exponents)
    return prod(factorial.(values(tmp)))
end

symmetric_polynomial_to_string(s::symmetric_polynomial) = is_elementary(s) ? "e"*string(sum(s.exponents)) : "S"*string(s.exponents)
function monomial_to_string(x::monomial)
    head = x.coeff>0 ? "  + " : "  - "
    head *= x.coeff.den == 1 ? string(abs(x.coeff.num)) : string(abs(x.coeff))
    for (s,e) in x.factors
        head *= " * "
        head *= e == 1 ? symmetric_polynomial_to_string(s) : symmetric_polynomial_to_string(s)*"^$e"
    end
    return head
end
polynomial_to_string(s::polynomial) = prod([monomial_to_string(x) for x in s.terms])

Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::polynomial) = join(io,polynomial_to_string(x))
