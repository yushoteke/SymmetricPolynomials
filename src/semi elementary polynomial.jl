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
    k.terms[x] = 1
    return k
end

function Base.push!(x::semi_elementary_polynomial{N}, k::semi_elementary_monomial{N}, v::Union{Integer,Rational}) where {N}
    if !haskey(x.terms, k)
        x.terms[k] = v
    elseif x.terms[k] != -v
        x.terms[k] += v
    else
        pop!(x.terms, k)
    end
end

Base.:(==)(x::semi_elementary_polynomial, y::semi_elementary_polynomial) = x.terms == y.terms

"""
    lower_order_representation(p,x)

Rewrite an admissible symmetric polynomial (as explained in the documentation for decompose)
In lower order terms, and add them to polynomial p.

For example, we want to write ``S(0,1,3)=x^3(y+z)+y^3(x+z)+z^3(x+y)`` into an equivalent representation,
but with terms of lower order. We would take an initial guess, ``S(0,1,3)\\approx S(0,1,2)S(0,0,1)``.

It turns out ``S(0,1,2)S(0,0,1)=S(0,1,3)+2S(0,2,2)+2S(1,1,1)(0,0,1)``, so if we rewrite it as
``S(0,1,3)=S(0,1,2)S(0,0,1)-2S(0,2,2)-2S(1,1,1)(0,0,1)`` every term has lower order.
"""
function lower_order_representation(polynomial::semi_elementary_polynomial{N}, x::semi_elementary_monomial{N}) where {N}
    is_elementary(x) && return
    sp = x.sp_term
    original_coefficient = polynomial.terms[x]
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

    f_keys, f_values = count_occurrences(factor)
    distribution_ways = ways_place_containers(f_values, num_highest)

    for j = 1:size(distribution_ways, 1)
        way = view(distribution_ways, j, :)
        representative = canonical_placement(f_values, way)
        new_term = TupleTools.sort(ntuple(i -> factor[i] + representative[i], N))
        coeff = 1
        for i = 2:length(f_keys)
            if f_keys[i-1] == f_keys[i] - 1
                coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        push!(polynomial, semi_elementary_monomial(new_term, x.powers), -coeff * original_coefficient)
    end
    tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], N)
    push!(polynomial, semi_elementary_monomial(factor, tmp), original_coefficient)
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
        lower_order_representation(polynomial, highest_term)
        highest_term = last(polynomial.terms).first
    end
    return polynomial
end

function to_string(x::semi_elementary_polynomial)
    head = ""
    for (i, j) in x.terms
        if j > 0
            head *= "+"
        end
        head *= string(j)
        head *= "*"
        head *= to_string(i)
    end
    return head
end
Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::semi_elementary_polynomial) = join(io, to_string(x))

function lower_order_representation2(polynomial::semi_elementary_polynomial{N},
    vec::Vector{semi_elementary_monomial{N}}) where {N}
    is_elementary(vec[1]) && return
    sp = vec[1].sp_term
    original_coefficients = map(x -> polynomial.terms[x], vec)
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

    f_keys, f_values = count_occurrences(factor)
    distribution_ways = ways_place_containers(f_values, num_highest)

    for j = 1:size(distribution_ways, 1)
        way = view(distribution_ways, j, :)
        representative = canonical_placement(f_values, way)
        new_term = TupleTools.sort(ntuple(i -> factor[i] + representative[i], N))
        coeff = 1
        for i = 2:length(f_keys)
            if f_keys[i-1] == f_keys[i] - 1
                coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        for (ind, x) ∈ enumerate(vec)
            push!(polynomial, semi_elementary_monomial(new_term, x.powers), -coeff * original_coefficients[ind])
        end
    end
    for (ind, x) ∈ enumerate(vec)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], N)
        push!(polynomial, semi_elementary_monomial(factor, tmp), original_coefficients[ind])
    end
end

function decompose2(x::semi_elementary_monomial{N}) where {N}
    p = semi_elementary_polynomial(N)
    p.terms[x] = 1
    while true
        is_elementary(last(p.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        highest_terms = [last(p.terms).first]
        cur_semi_token = lastindex(p.terms)
        while cur_semi_token != startof(p.terms)
            cur_semi_token = regress((p.terms, cur_semi_token))
            key = deref_key((p.terms, cur_semi_token))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end
        lower_order_representation2(p, highest_terms)

    end
    return p
end

function decompose3(x::semi_elementary_monomial{N}) where {N}
    p = semi_elementary_polynomial(N)
    p.terms[x] = 1
    while true
        is_elementary(last(p.terms).first) && break

        #combines decompose2 and decompose strategies
        st = lastindex(p.terms)
        if st == startof(p.terms)
            lower_order_representation(p, deref_key((p.terms, st)))
            continue
        end

        prev_st = regress((p.terms, st))
        if deref_key((p.terms, st)).sp_term != deref_key((p.terms, prev_st)).sp_term
            lower_order_representation(p, deref_key((p.terms, st)))
            continue
        end

        highest_terms = [last(p.terms).first]
        while st != startof(p.terms)
            st = regress((p.terms, st))
            key = deref_key((p.terms, st))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end
        lower_order_representation2(p, highest_terms)

    end
    return p
end

function lower_order_representation3(polynomial::semi_elementary_polynomial{N},
    vec::Vector{semi_elementary_monomial{N}}) where {N}
    is_elementary(vec[1]) && return
    sp = vec[1].sp_term
    original_coefficients = map(x -> polynomial.terms[x], vec)
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

    f_keys, f_values = count_occurrences(factor)
    distribution_ways_iter = ways_place_containers_iterator(f_values, num_highest)

    for way ∈ distribution_ways_iter
        representative = canonical_placement(f_values, way)
        new_term = TupleTools.sort(ntuple(i -> factor[i] + representative[i], N))
        coeff = 1
        for i = 2:length(f_keys)
            if f_keys[i-1] == f_keys[i] - 1
                coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        for (ind, x) ∈ enumerate(vec)
            push!(polynomial, semi_elementary_monomial(new_term, x.powers), -coeff * original_coefficients[ind])
        end
    end
    for (ind, x) ∈ enumerate(vec)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], N)
        push!(polynomial, semi_elementary_monomial(factor, tmp), original_coefficients[ind])
    end
end

function decompose4(x::semi_elementary_monomial{N}) where {N}
    p = semi_elementary_polynomial(N)
    p.terms[x] = 1
    while true
        is_elementary(last(p.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        highest_terms = [last(p.terms).first]
        cur_semi_token = lastindex(p.terms)
        while cur_semi_token != startof(p.terms)
            cur_semi_token = regress((p.terms, cur_semi_token))
            key = deref_key((p.terms, cur_semi_token))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end
        lower_order_representation3(p, highest_terms)

    end
    return p
end

function decompose5(x::semi_elementary_monomial{N}) where {N}
    p = semi_elementary_polynomial(N)
    p.terms[x] = 1
    highest_terms = semi_elementary_monomial{N}[]
    while true
        is_elementary(last(p.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        empty!(highest_terms)
        push!(highest_terms, last(p.terms).first)
        cur_semi_token = lastindex(p.terms)
        while cur_semi_token != startof(p.terms)
            cur_semi_token = regress((p.terms, cur_semi_token))
            key = deref_key((p.terms, cur_semi_token))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end
        lower_order_representation3(p, highest_terms)

    end
    return p
end

function lower_order_representation4(polynomial::semi_elementary_polynomial{N},
    vec::Vector{semi_elementary_monomial{N}},
    iter::ways_place_containers_iterator) where {N}
    is_elementary(vec[1]) && return
    sp = vec[1].sp_term
    original_coefficients = map(x -> polynomial.terms[x], vec)
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

    f_keys, f_values = count_occurrences(factor)
    reset!(iter, f_values, num_highest)

    for way ∈ iter
        representative = canonical_placement(f_values, way)
        new_term = TupleTools.sort(ntuple(i -> factor[i] + representative[i], N))
        coeff = 1
        for i = 2:length(f_keys)
            if f_keys[i-1] == f_keys[i] - 1
                coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        for (ind, x) ∈ enumerate(vec)
            push!(polynomial, semi_elementary_monomial(new_term, x.powers), -coeff * original_coefficients[ind])
        end
    end
    for (ind, x) ∈ enumerate(vec)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], N)
        push!(polynomial, semi_elementary_monomial(factor, tmp), original_coefficients[ind])
    end
end

function decompose6(x::semi_elementary_monomial{N}) where {N}
    p = semi_elementary_polynomial(N)
    p.terms[x] = 1
    iter = ways_place_containers_iterator([1], 1)
    highest_terms = semi_elementary_monomial{N}[]
    while true
        is_elementary(last(p.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        empty!(highest_terms)
        push!(highest_terms, last(p.terms).first)
        cur_semi_token = lastindex(p.terms)
        while cur_semi_token != startof(p.terms)
            cur_semi_token = regress((p.terms, cur_semi_token))
            key = deref_key((p.terms, cur_semi_token))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end
        lower_order_representation4(p, highest_terms, iter)

    end
    return p
end

"""
    Structure to hold all the arrays that are used in
    intermediate computation to avoid reallocating

"""
struct polynomial_decomposer{N}
    iter::ways_place_containers_iterator
    f_keys::Vector{Int64}
    f_values::Vector{Int64}
    representative::Vector{Int64}
    new_term::Vector{Int64}
    factor::Vector{Int64}
    original_coefficients::Vector{Int128}
    highest_terms::Vector{semi_elementary_monomial{N}}
    polynomial::semi_elementary_polynomial{N}
end

polynomial_decomposer(N) = polynomial_decomposer{N}(ways_place_containers_iterator([1], 1), [], [], [], [], [], [], [], semi_elementary_polynomial(N))

function decompose7(x::semi_elementary_monomial{N}) where {N}
    D = polynomial_decomposer(N)
    D.polynomial.terms[x] = 1
    while true
        is_elementary(last(D.polynomial.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        empty!(D.highest_terms)
        push!(D.highest_terms, last(D.polynomial.terms).first)
        cur_semi_token = lastindex(D.polynomial.terms)
        while cur_semi_token != startof(D.polynomial.terms)
            cur_semi_token = regress((D.polynomial.terms, cur_semi_token))
            key = deref_key((D.polynomial.terms, cur_semi_token))
            key.sp_term == D.highest_terms[1].sp_term ? push!(D.highest_terms, key) : break
        end
        lower_order_representation5(D)

    end
    return D.polynomial
end

function lower_order_representation5(D::polynomial_decomposer{N}) where {N}
    is_elementary(D.highest_terms[1]) && return
    sp = D.highest_terms[1].sp_term
    empty!(D.original_coefficients)
    for x in D.highest_terms
        push!(D.original_coefficients, D.polynomial.terms[x])
    end
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    empty!(D.factor)
    for i = 1:N
        push!(D.factor, i > N - num_highest ? sp[i] - 1 : sp[i])
    end

    count_occurrences!(D.factor, D.f_keys, D.f_values)
    reset!(D.iter, D.f_values, num_highest)

    for way ∈ D.iter
        canonical_placement!(D.f_values, way, D.representative)
        empty!(D.new_term)
        for i = 1:N
            push!(D.new_term, D.factor[i] + D.representative[i])
        end
        coeff = 1
        for i = 2:length(D.f_keys)
            if D.f_keys[i-1] == D.f_keys[i] - 1
                coeff *= binomial(D.f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        for (ind, x) ∈ enumerate(D.highest_terms)
            tup = ntuple(i -> D.new_term[i], Val(N))
            push!(D.polynomial, semi_elementary_monomial(tup, x.powers), -coeff * D.original_coefficients[ind])
        end
    end
    for (ind, x) ∈ enumerate(D.highest_terms)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
        factor = ntuple(i -> D.factor[i], Val(N))
        push!(D.polynomial, semi_elementary_monomial(factor, tmp), D.original_coefficients[ind])
    end
end

function decompose8(x::semi_elementary_monomial{N}) where {N}
    D = polynomial_decomposer(N)
    D.polynomial.terms[x] = 1
    while true
        is_elementary(last(D.polynomial.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        empty!(D.highest_terms)
        push!(D.highest_terms, last(D.polynomial.terms).first)
        cur_semi_token = lastindex(D.polynomial.terms)
        while cur_semi_token != startof(D.polynomial.terms)
            cur_semi_token = regress((D.polynomial.terms, cur_semi_token))
            key = deref_key((D.polynomial.terms, cur_semi_token))
            key.sp_term == D.highest_terms[1].sp_term ? push!(D.highest_terms, key) : break
        end
        lower_order_representation6(D)

    end
    return D.polynomial
end

function lower_order_representation6(D::polynomial_decomposer{N}) where {N}
    is_elementary(D.highest_terms[1]) && return
    sp = D.highest_terms[1].sp_term
    empty!(D.original_coefficients)
    for x in D.highest_terms
        push!(D.original_coefficients, D.polynomial.terms[x])
    end
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    empty!(D.factor)
    for i = 1:N
        push!(D.factor, i > N - num_highest ? sp[i] - 1 : sp[i])
    end

    count_occurrences!(D.factor, D.f_keys, D.f_values)
    reset!(D.iter, D.f_values, num_highest)

    for way ∈ D.iter
        canonical_placement!(D.f_values, way, D.representative)
        empty!(D.new_term)
        for i = 1:N
            push!(D.new_term, D.factor[i] + D.representative[i])
        end
        coeff = 1
        for i = 2:length(D.f_keys)
            if D.f_keys[i-1] == D.f_keys[i] - 1
                coeff *= binomial(D.f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        for (ind, x) ∈ enumerate(D.highest_terms)
            tup = ntuple(i -> D.new_term[i], Val(N))
            push!(D.polynomial, semi_elementary_monomial(tup, x.powers, false), -coeff * D.original_coefficients[ind])
        end
    end
    for (ind, x) ∈ enumerate(D.highest_terms)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
        factor = ntuple(i -> D.factor[i], Val(N))
        push!(D.polynomial, semi_elementary_monomial(factor, tmp, false), D.original_coefficients[ind])
    end
end

function decompose9(x::semi_elementary_monomial{N}) where {N}
    D = polynomial_decomposer(N)
    D.polynomial.terms[x] = 1
    while true
        is_elementary(last(D.polynomial.terms).first) && break

        #find all semi elementary monomials whose non elementary part are same
        empty!(D.highest_terms)
        push!(D.highest_terms, last(D.polynomial.terms).first)
        cur_semi_token = lastindex(D.polynomial.terms)
        while cur_semi_token != startof(D.polynomial.terms)
            cur_semi_token = regress((D.polynomial.terms, cur_semi_token))
            key = deref_key((D.polynomial.terms, cur_semi_token))
            @inbounds key.sp_term == D.highest_terms[1].sp_term ? push!(D.highest_terms, key) : break
        end
        lower_order_representation7(D)

    end
    return D.polynomial
end

function lower_order_representation7(D::polynomial_decomposer{N}) where {N}
    @inbounds is_elementary(D.highest_terms[1]) && return
    @inbounds sp = D.highest_terms[1].sp_term
    empty!(D.original_coefficients)
    for x in D.highest_terms
        @inbounds push!(D.original_coefficients, D.polynomial.terms[x])
    end
    num_highest = N - findfirst(i -> i == sp[end], sp) + 1
    empty!(D.factor)
    for i = 1:N
        @inbounds push!(D.factor, i > N - num_highest ? sp[i] - 1 : sp[i])
    end

    count_occurrences!(D.factor, D.f_keys, D.f_values)
    reset!(D.iter, D.f_values, num_highest)

    for way ∈ D.iter
        canonical_placement!(D.f_values, way, D.representative)
        empty!(D.new_term)
        for i = 1:N
            push!(D.new_term, D.factor[i] + D.representative[i])
        end
        coeff = 1
        for i = 2:length(D.f_keys)
            if D.f_keys[i-1] == D.f_keys[i] - 1
                @inbounds coeff *= binomial(D.f_values[i] - way[i] + way[i-1], way[i-1])
            end
        end
        sort!(D.new_term)
        @inbounds tup = ntuple(i -> D.new_term[i], Val(N))
        for (ind, x) ∈ enumerate(D.highest_terms)
            @inbounds push!(D.polynomial, semi_elementary_monomial(tup, x.powers, true), -coeff * D.original_coefficients[ind])
        end
    end
    sort!(D.factor)
    @inbounds tup = ntuple(i -> D.factor[i], Val(N))
    for (ind, x) ∈ enumerate(D.highest_terms)
        @inbounds tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
        @inbounds push!(D.polynomial, semi_elementary_monomial(tup, tmp, true), D.original_coefficients[ind])
    end
end