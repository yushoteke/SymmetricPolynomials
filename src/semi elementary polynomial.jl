struct semi_elementary_polynomial{N}
    terms::SortedDict{semi_elementary_monomial{N},Int64,Base.Order.ForwardOrdering}
    semi_elementary_polynomial(N) = new{N}(SortedDict{semi_elementary_monomial{N},Int64}())
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

function lower_order_representation(polynomial::semi_elementary_polynomial{N},x::semi_elementary_monomial{N}) where {N}
    is_elementary(x) && return
    sp = x.sp_term
    original_coefficient = polynomial.terms[x]
    Z = ntuple(i->0,N)
    num_highest = N - findfirst(i->i==sp[end],sp) + 1
    factor1 = ntuple(i->i > N - num_highest ? 1 : 0,N)
    factor2 = sp .- factor1

    f2_containers = counter(factor2)
    f2_keys = sort(collect(keys(f2_containers)))
    f2_values = map(i->f2_containers[i],f2_keys)
    distribution_ways = ways_place_containers(f2_values,num_highest)

    for way in distribution_ways
        representative = canonical_placement(f2_values,way)
        new_term = sorttuple(add_tuple_array(factor2,representative))
        coeff = 1
        for i=2:length(f2_keys)
            if f2_keys[i-1]==f2_keys[i]-1
                coeff *= binomial(f2_values[i]-way[i]+way[i-1],way[i-1])
            end
        end
        #=
        if new_term[end]!=1
            new_monomial = semi_elementary_monomial(new_term,x.powers)
        else
            tmp1::NTuple{N,eltype(factor2)} = elementary_representation(new_term)
            new_monomial = semi_elementary_monomial(Z,tmp1.+x.powers)
        end
        =#
        new_monomial::semi_elementary_monomial{N} = new_term[end]!=1 ?
                semi_elementary_monomial(new_term,x.powers) : semi_elementary_monomial(Z,elementary_representation(new_term).+x.powers)
                                                
        push!(polynomial,new_monomial,-coeff * original_coefficient)
    end
    if factor2[end]==1
        tmp::NTuple{N,eltype(factor2)} = elementary_representation(factor1) .+ elementary_representation(factor2) .+ x.powers
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
        if j>0
            head *= "+"
        end
        head *= string(j)
        head *= "*"
        head *= to_string(i)
    end
    return head
end
Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::semi_elementary_polynomial) = join(io,to_string(x))
