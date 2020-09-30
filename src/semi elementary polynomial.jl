"""
    semi_elementary_polynomial{N}

A polynomial with integer coefficient, where each term is a semi elementary monomial.
Internally it's represented as a sorted dict, which allows fast lookup of the term
with the highest order, and fast insertion.
"""
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

"""
    lower_order_representation(p,x)

Rewrite an admissible symmetric polynomial (as explained in the documentation for decompose)
In lower order terms, and add them to polynomial p.

For example, we want to write ``S(0,1,3)=x^3(y+z)+y^3(x+z)+z^3(x+y)`` into an equivalent representation,
but with terms of lower order. We would take an initial guess, ``S(0,1,3)\\approx S(0,1,2)S(0,0,1)``.

It turns out ``S(0,1,2)S(0,0,1)=S(0,1,3)+2S(0,2,2)+2S(1,1,1)(0,0,1)``, so if we rewrite it as
``S(0,1,3)=S(0,1,2)S(0,0,1)-2S(0,2,2)-2S(1,1,1)(0,0,1)`` every term has lower order.
"""
function lower_order_representation(polynomial::semi_elementary_polynomial{N},x::semi_elementary_monomial{N}) where {N}
    is_elementary(x) && return
    sp = x.sp_term
    original_coefficient = polynomial.terms[x]
    num_highest = N - findfirst(i->i==sp[end],sp) + 1
    factor = ntuple(i->i > N - num_highest ? sp[i] - 1 : sp[i],N)

    f_keys,f_values = count_occurrences(factor)
    distribution_ways = ways_place_contain
push!(LOAD_PATH,"../src/")ers(f_values,num_highest)

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

"""
    decompose(x)

Rewrites a symmetric polynomial as a polynomial of elementary symmetric polynomials.
This function only factors symmetric polynomials where every term has the same order
and the powers has the same partition.
For example, ``x^3(y+z)+y^3(x+z)+z^3(x+y)`` is ok. Because the power of every term
is ``(0,1,3)``.

``(x^2+y^2+z^2)+(x^3+y^3+z^3)`` on the other hand is not ok, because some terms' power
are ``(0,0,2)`` while others are ``(0,0,3)``.
"""
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
