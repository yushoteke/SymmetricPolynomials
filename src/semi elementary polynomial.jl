struct semi_elementary_polynomial{N}
    terms::SortedDict{semi_elementary_monomial{N},Int128,Base.Order.ForwardOrdering}
    semi_elementary_polynomial(N) = new{N}(SortedDict{semi_elementary_monomial{N},Int128}())
end

function to_polynomial(x::semi_elementary_monomial{N}) where {N}
    k = semi_elementary_polynomial(N)
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
    num_highest = N - findfirst(i->i==sp[end],sp) + 1
    factor = ntuple(i->i > N - num_highest ? sp[i] - 1 : sp[i],N)

    f_keys,f_values = count_occurrences(factor)
    distribution_ways = ways_place_containers(f_values,num_highest)

    for j=1:size(distribution_ways,1)
        way = view(distribution_ways,j,:)
        representative = canonical_placement(f_values,way)
        new_term = TupleTools.sort(ntuple(i->factor[i] + representative[i],N))
        coeff = 1
        for i=2:length(f_keys)
            if f_keys[i-1]==f_keys[i]-1
                coeff *= binomial(f_values[i] - way[i] + way[i-1],way[i-1])
            end
        end
        push!(polynomial,semi_elementary_monomial(new_term,x.powers),-coeff * original_coefficient)
    end
    tmp = ntuple(i->i==num_highest ? x.powers[i] + 1 : x.powers[i],N)
    push!(polynomial,semi_elementary_monomial(factor,tmp),original_coefficient)
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
