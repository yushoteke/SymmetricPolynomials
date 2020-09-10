#=
function ordering_m1_n0(m,n)
    #generates all unique orderings or m 1s and n 0s
    #=
    For Example, 2 1s and 2 0s
    1100,1010,1001,0110,0101,0011
    =#
    #base case
    if m<=n<=0
        return Array{Int64,1}[]
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
=#
function ordering_m1_n0(m,n)
    arr = zeros(Int,binomial(m+n,m),m+n)
    set_1s_0s(arr,1,1,m+n,m,n)
    return arr
end
function set_1s_0s(arr,h,w,wmax,m,n)
    x = binomial(m+n-1,n)
    arr[h:h+x-1,w] .= 1
    if w!=wmax
        set_1s_0s(arr,h,w+1,wmax,m-1,n)
        set_1s_0s(arr,h+x,w+1,wmax,m,n-1)
    end
end

function ways_place_containers(x,m)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with
    #positive capacity x_i. Now we want to place m balls into these containers
    #without breaking the capacity

    #base case, assume x <= sum(x), so there is at least 1 valid solution
    length(x)==1 && return Array{Int,1}[[m]]
    results = Array{Int,1}[]
    maximum = min(x[end],m)
    minimum = max(0,m - sum(x[1:end-1]))
    for balls_in_last_container = minimum:maximum
        subsolutions = ways_place_containers(x[1:end-1],m-balls_in_last_container)
        if subsolutions!= nothing
            for i in subsolutions
                push!(i,balls_in_last_container)
            end
        end
        append!(results,subsolutions)
    end
    return results
end

function canonical_placement(x,y)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with positive capacity x_i.
    #assume y=[y_1,y_2,y_3...y_n] means placed 0<=y_i<=x_i balls into container i
    #returns one way this could have happened.
    #For example, x=[1,2,6],y=[1,1,2] could be 1|10|010100 or 1|01|110000
    tmp = zeros(Int,sum(x))
    starting_locations = cumsum([1,x...])
    for (i,x) in enumerate(starting_locations[1:end-1])
        tmp[x:x-1+y[i]] .= 1
    end
    return tmp
end

add_tuple_array(x::NTuple{N,T},y) where {N,T<:Any} = ntuple(i->x[i]+y[i],N)


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

elementary_symmetric_polynomial_to_string(x::elementary_symmetric_polynomial) = "e"*string(x.order)
function elementary_monomial_to_string(x::elementary_monomial)
    head = x.coeff>0 ? "  +" : "  -"
    head *= x.coeff.den == 1 ? string(abs(x.coeff.num)) : string(abs(x.coeff))
    head *= "*"
    for (i,n) in enumerate(x.exponents)
        if n == 0
            continue
        elseif n == 1
            head *= "e$i "
        else
            head *= "e$i^$n "
        end
    end
    return head
end
elementary_polynomial_to_string(x::elementary_polynomial) = prod([elementary_monomial_to_string(i) for i in x.terms])

Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::elementary_polynomial) = join(io,elementary_polynomial_to_string(x))
#=
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
polynomial_to_string(x::polynomial) = prod([monomial_to_string(i) for i in x.terms])

Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::polynomial) = join(io,polynomial_to_string(x))
=#
