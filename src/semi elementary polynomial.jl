
"""
A semi elementary polynomial is a polynomial with semi elementary monomial terms
It is represented as a SortedDict because this allows
1.fast merging terms
2.Fast peek (finding the highest order term to lower first)
"""
struct semi_elementary_polynomial{N}
    terms::SortedDict{semi_elementary_monomial{N},Rational{BigInt},Base.Order.ForwardOrdering}
    semi_elementary_polynomial(N) = new{N}(SortedDict{semi_elementary_monomial{N},Rational{BigInt}}())

end

function to_polynomial(x::semi_elementary_monomial)
    k = semi_elementary_polynomial(dim(x))
    k.terms[x]=1
    return k
end

function Base.push!(x::semi_elementary_polynomial{N},k::semi_elementary_monomial{N},v::Union{Integer,Rational}) where {N}
    if !haskey(x.terms,k)
        x.terms[k] = v
    elseif x.terms[k] != -v
        x.terms[k] += v
    else
        pop!(x.terms,k)
    end
end

Base.:(==)(x::semi_elementary_polynomial,y::semi_elementary_polynomial) = x.terms == y.terms

function muladd!(x::semi_elementary_polynomial{N},y::semi_elementary_polynomial{N},c::Union{Integer,Rational}) where {N}
    for (i,j) in y.terms
        push!(x,i,c*j)
    end
end

function lower_order_representation(polynomial::semi_elementary_polynomial{N},x::semi_elementary_monomial{N}) where {N}
    #pick the key x from polynomial
    #replace that term with an equivalent representation, but of lower order,
    #and store it into polynomial in place
    is_elementary(x) && return
    sp = x.sp_term
    original_coefficient = polynomial.terms[x]
    Z = ntuple(i->0,N)
    num_highest = N - findfirst(i->i==sp[end],sp) + 1
    factor1 = ntuple(i->i > N - num_highest ? 1 : 0,N)
    factor2 = sp .- factor1

    c = symmetry_factor(factor2)

    #e.g S(3,3,3,3,3,3,4,4,5) -> [6,2,1]
    container_sizes = count_occurrences(sp)
    #eg [6,2,1],4 -> [4,0,0],[3,1,0],[3,0,1]...
    distribution_ways = ways_place_containers(container_sizes,num_highest)
    for way in distribution_ways
        representative = canonical_placement(container_sizes,way)

        new_term = sorttuple(add_tuple_array(factor2,representative))
        coeff = symmetry_factor(new_term) * prod([binomial(container_sizes[i],way[i]) for i=1:length(way)]) // c

        if new_term[end]!=1
            new_monomial = semi_elementary_monomial(new_term,x.powers)
        else
            new_monomial = semi_elementary_monomial(Z,elementary_representation(new_term).+x.powers)
        end
        push!(polynomial,new_monomial,-coeff * original_coefficient)
        #add a negative sign here because later we will move almost everything to the other side of the equation
    end
    #turn the product of factor1 and factor2 into a monomial
    if factor2[end]==1
        tmp = elementary_representation(factor1) .+ elementary_representation(factor2) .+ x.powers
        push!(polynomial,semi_elementary_monomial(Z,tmp),original_coefficient)
    else
        push!(polynomial,semi_elementary_monomial(factor2,elementary_representation(factor1).+x.powers),original_coefficient)
    end
end

function decompose(x::semi_elementary_monomial{N}) where {N}
    polynomial = semi_elementary_polynomial(N)
    polynomial.terms[x] = 1
    highest_term = last(polynomial.terms).first
    while !is_elementary(highest_term)
        #println("highest term is "*string(highest_term))
        lower_order_representation(polynomial,highest_term)
        highest_term = last(polynomial.terms).first
    end
    return polynomial
end


function to_string(x::semi_elementary_polynomial)
    head = ""
    for (i,j) in x.terms
        head *= j>0 ? " +" : " -"
        head *= j.den==1 ? string(abs(j.num)) : string(abs(j))
        head *= "*"
        head *= to_string(i)
    end
    return head
end
Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::semi_elementary_polynomial) = join(io,to_string(x))
