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
            key.sp_term == D.highest_terms[1].sp_term ? push!(D.highest_terms, key) : break
        end
        lower_order_representation7(D)

    end
    return D.polynomial
end

function lower_order_representation7(D::polynomial_decomposer{N}) where {N}
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
        sort!(D.new_term)
        tup = ntuple(i -> D.new_term[i], Val(N))
        for (ind, x) ∈ enumerate(D.highest_terms)
            push!(D.polynomial, semi_elementary_monomial(tup, x.powers, true), -coeff * D.original_coefficients[ind])
        end
    end
    sort!(D.factor)
    tup = ntuple(i -> D.factor[i], Val(N))
    for (ind, x) ∈ enumerate(D.highest_terms)
        tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
        push!(D.polynomial, semi_elementary_monomial(tup, tmp, true), D.original_coefficients[ind])
    end
end

struct semi_elementary_polynomial2{N}
    terms::Dict{semi_elementary_monomial{N},Int128}
    sort_terms::SortedSet{semi_elementary_monomial{N}}
    semi_elementary_polynomial2(N) = new{N}(Dict{semi_elementary_monomial{N},Int128}(),
        SortedSet{semi_elementary_monomial{N}}())
end

function Base.push!(x::semi_elementary_polynomial2{N}, k::semi_elementary_monomial{N}, v::Union{Integer,Rational}) where {N}
    if !haskey(x.terms, k)
        x.terms[k] = v
        Base.insert!(x.sort_terms, k)
    elseif x.terms[k] != -v
        x.terms[k] += v
    else
        pop!(x.terms, k)
        pop!(x.sort_terms, k)
    end
end

function poly2_to_poly(x::semi_elementary_polynomial2{N}) where {N}
    p = semi_elementary_polynomial(N)
    for (k, v) ∈ x.terms
        p.terms[k] = v
    end
    return p
end

function decompose10(x::semi_elementary_monomial{N}) where {N}
    polynomial = semi_elementary_polynomial2(N)
    push!(polynomial, x, 1)

    iter = ways_place_containers_iterator([1], 1) #init variables are arbitrary, since they will be replaced anyways
    highest_terms = semi_elementary_monomial{N}[]
    original_coefficients = Int128[]
    f_keys = Int64[]
    f_values = Int64[]
    representative = Int64[]

    while true
        is_elementary(last(polynomial.sort_terms)) && break

        #setup part, find highest terms
        empty!(highest_terms)
        push!(highest_terms, last(polynomial.sort_terms))
        cur_semi_token = lastindex(polynomial.sort_terms)
        while cur_semi_token != startof(polynomial.sort_terms)
            cur_semi_token = regress((polynomial.sort_terms, cur_semi_token))
            key = deref((polynomial.sort_terms, cur_semi_token))
            key.sp_term == highest_terms[1].sp_term ? push!(highest_terms, key) : break
        end


        #lower polynomial order part
        sp = highest_terms[1].sp_term
        empty!(original_coefficients)
        for term in highest_terms
            push!(original_coefficients, polynomial.terms[term])
        end
        num_highest = N - findfirst(i -> i == sp[end], sp) + 1
        factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

        count_occurrences!(factor, f_keys, f_values)
        reset!(iter, f_values, num_highest)

        for way ∈ iter
            canonical_placement!(f_values, way, representative)
            new_term = ntuple(i -> factor[i] + representative[i], Val(N))
            coeff = 1
            for i = 2:length(f_keys)
                if f_keys[i-1] == f_keys[i] - 1
                    coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                end
            end
            for (ind, x) ∈ enumerate(highest_terms)
                push!(polynomial, semi_elementary_monomial(new_term, x.powers), -coeff * original_coefficients[ind])
            end
        end

        for (ind, x) ∈ enumerate(highest_terms)
            tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
            push!(polynomial, semi_elementary_monomial(factor, tmp), original_coefficients[ind])
        end
    end
    return polynomial
end

#for decompose 11, set 2 polynomials, where one contains only elementary elements

#since sorting is a big issue, we can create many small sorted containers of different "rank".
#for example, suppose we want to decompose S((0,0,0,0,0,80),(0,0,0,0,0,0))
#Then we could create a length 80 vector V, where V[i] can only hold elements of the form 
#S((0,x,x,x,x,i),(x,x,x,x,x,x)). This way, the size of each sorted container is considerably smaller
#this method also incorporates the previous idea

#try @inline to improve code readability
using Combinatorics

struct semi_elementary_polynomial3{N}
    coeffs::Dict{semi_elementary_monomial{N},Int128}
    s_terms::Dict{NTuple{N,Int64},Set{semi_elementary_monomial{N}}}
end

function semi_elementary_polynomial3(N)
    return semi_elementary_polynomial3{N}(Dict(), Dict())
end

function poly3_to_poly(x::semi_elementary_polynomial3{N}) where {N}
    p = semi_elementary_polynomial(N)
    for (k, v) ∈ x.coeffs
        p.terms[k] = v
    end
    return p
end

function Base.push!(p::semi_elementary_polynomial3{N}, x::semi_elementary_monomial{N}, v::Union{Integer,Rational}) where {N}
    if !haskey(p.coeffs, x)
        p.coeffs[x] = v
        if !haskey(p.s_terms, x.sp_term)
            p.s_terms[x.sp_term] = Set()
        end
        if x.sp_term[N] != 0
            push!(p.s_terms[x.sp_term], x)
        end
    elseif p.coeffs[x] != -v
        p.coeffs[x] += v
    else
        pop!(p.coeffs, x)
        if x.sp_term[N] != 0
            pop!(p.s_terms[x.sp_term], x)
        end
    end
end

function push_timeit!(p::semi_elementary_polynomial3{N}, x::semi_elementary_monomial{N}, v::Union{Integer,Rational}, timer) where {N}


    if !haskey(p.coeffs, x)
        @timeit timer "set coeff" p.coeffs[x] = v
        if !haskey(p.s_terms, x.sp_term)
            @timeit timer "create new set" p.s_terms[x.sp_term] = Set()
            p.calls[1] += 1
        end
        if x.sp_term[N] != 0
            tmp = @timeit timer "get container alias" p.s_terms[x.sp_term]
            @timeit timer "push new term" push!(tmp, x)
        end
    elseif p.coeffs[x] != -v
        @timeit timer "set coeff" p.coeffs[x] += v
    else
        @timeit timer "pop coeff" pop!(p.coeffs, x)
        if x.sp_term[N] != 0
            @timeit timer "pop term" pop!(p.s_terms[x.sp_term], x)
        end
    end

end

function partition_to_tuple(p::Vector{Int64}, ::Val{N}) where {N}
    # for example, [4,2,1], 5 returns (0,0,1,2,4)
    #which means 
    #1 -> 0
    #2 -> 0
    #3 -> p[3]
    #4 -> p[2]
    #5 -> p[1]
    return ntuple(i -> (i >= N + 1 - length(p)) ? p[N+1-i] : 0, N)
end

function decompose11(x::semi_elementary_monomial{N}) where {N}
    rank = sum(x.sp_term)
    p = semi_elementary_polynomial3(N)
    push!(p, x, 1)

    iter = ways_place_containers_iterator([1], 1)
    f_keys = Int64[]
    f_values = Int64[]
    representative = Int64[]
    sp_iter = sp_term_iterator(rank, Val(N))

    for sp ∈ sp_iter
        highest_terms = p.s_terms[sp] |> collect
        orig_coeffs = map(x -> p.coeffs[x], highest_terms)
        num_highest = N - findfirst(i -> i == sp[end], sp) + 1
        factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

        count_occurrences!(factor, f_keys, f_values)
        reset!(iter, f_values, num_highest)

        for way ∈ iter
            canonical_placement!(f_values, way, representative)
            new_term = ntuple(i -> factor[i] + representative[i], N)
            coeff = 1
            for i = 2:length(f_keys)
                if f_keys[i-1] == f_keys[i] - 1
                    coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                end
            end
            for (ind, x) ∈ enumerate(highest_terms)
                push!(p, semi_elementary_monomial(new_term, x.powers), -coeff * orig_coeffs[ind])
            end
        end

        for (ind, x) ∈ enumerate(highest_terms)
            tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
            push!(p, semi_elementary_monomial(factor, tmp), orig_coeffs[ind])
        end
        #pop!(p.s_terms, sp)
    end
    return p
end

function decompose12(x::semi_elementary_monomial{N}) where {N}
    rank = sum(x.sp_term)
    p = semi_elementary_polynomial3(N)
    push!(p, x, 1)

    iter = ways_place_containers_iterator([1], 1)
    f_keys = Int64[]
    f_values = Int64[]
    representative = Int64[]
    sp_iter = sp_term_iterator(rank, Val(N))

    for sp ∈ sp_iter
        highest_terms = p.s_terms[sp] |> collect
        orig_coeffs = map(x -> p.coeffs[x], highest_terms)
        num_highest = N - findfirst(i -> i == sp[end], sp) + 1
        factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

        count_occurrences!(factor, f_keys, f_values)
        reset!(iter, f_values, num_highest)

        for way ∈ iter
            canonical_placement!(f_values, way, representative)
            new_term = ntuple(i -> factor[i] + representative[i], N)
            coeff = 1
            for i = 2:length(f_keys)
                if f_keys[i-1] == f_keys[i] - 1
                    coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                end
            end
            for (ind, x) ∈ enumerate(highest_terms)
                push!(p, semi_elementary_monomial(new_term, x.powers), -coeff * orig_coeffs[ind])
            end
        end

        for (ind, x) ∈ enumerate(highest_terms)
            tmp = ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
            push!(p, semi_elementary_monomial(factor, tmp), orig_coeffs[ind])
        end
        pop!(p.s_terms, sp)
    end
    return p
end

using TimerOutputs

function decompose12_timeit(x::semi_elementary_monomial{N}) where {N}
    timer = TimerOutput()
    rank = sum(x.sp_term)
    p = semi_elementary_polynomial3(N)
    sizehint!(p.coeffs, 2^18)
    sizehint!(p.s_terms, 2^16)
    push!(p, x, 1)

    iter = ways_place_containers_iterator([1], 1)
    f_keys = Int64[]
    f_values = Int64[]
    representative = Int64[]
    sp_iter = sp_term_iterator(rank, Val(N))

    for sp ∈ sp_iter
        @timeit timer "collect" begin
            highest_terms = @timeit timer "collect highest term" p.s_terms[sp] |> collect
            orig_coeffs = @timeit timer "collect coefficients" map(x -> p.coeffs[x], highest_terms)
            num_highest = N - findfirst(i -> i == sp[end], sp) + 1
            factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)
        end

        count_occurrences!(factor, f_keys, f_values)
        reset!(iter, f_values, num_highest)

        @timeit timer "lower order" begin
            for way ∈ iter
                canonical_placement!(f_values, way, representative)

                new_term = @inbounds @timeit timer "create tuple" ntuple(i -> factor[i] + representative[i], Val(N))
                coeff = 1
                for i = 2:length(f_keys)
                    if f_keys[i-1] == f_keys[i] - 1
                        coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                    end
                end
                @timeit timer "push1" begin
                    for (ind, x) ∈ enumerate(highest_terms)
                        push_timeit!(p, semi_elementary_monomial(new_term, x.powers), -coeff * orig_coeffs[ind], timer)
                    end
                end
            end

            @timeit timer "push2" begin
                for (ind, x) ∈ enumerate(highest_terms)
                    tmp = @timeit timer "create tuple 2" ntuple(i -> i == num_highest ? x.powers[i] + 1 : x.powers[i], Val(N))
                    push_timeit!(p, semi_elementary_monomial(factor, tmp), orig_coeffs[ind], timer)
                end
            end
        end
        pop!(p.s_terms, sp)
    end
    println(p.calls[1])
    show(timer)
    return p
end

struct semi_elementary_polynomial4{N}
    coeffs::DefaultDict{semi_elementary_monomial{N},Int128,Int128}
    s_terms::Dict{NTuple{N,Int64},Vector{semi_elementary_monomial{N}}}
end

function semi_elementary_polynomial4(N)
    return semi_elementary_polynomial4{N}(DefaultDict{semi_elementary_monomial{N},Int128,Int128}(zero(Int128)), Dict())
end

function poly4_to_poly(x::semi_elementary_polynomial4{N}) where {N}
    p = semi_elementary_polynomial(N)
    for (k, v) ∈ x.coeffs
        p.terms[k] = v
    end
    return p
end

function Base.push!(p::semi_elementary_polynomial4{N}, x::semi_elementary_monomial{N}, v::Union{Integer,Rational}) where {N}
    if !haskey(p.coeffs, x)
        p.coeffs[x] = v
        !haskey(p.s_terms, x.sp_term) && (p.s_terms[x.sp_term] = [])
        x.sp_term[N] != 0 && push!(p.s_terms[x.sp_term], x)
    elseif p.coeffs[x] != -v
        p.coeffs[x] += v
    else
        pop!(p.coeffs, x)
    end
end

function decompose13(x::semi_elementary_monomial{N}) where {N}
    rank = sum(x.sp_term)
    p = semi_elementary_polynomial4(N)
    push!(p, x, 1)

    iter = ways_place_containers_iterator([1], 1)
    f_keys = Int64[]
    f_values = Int64[]
    representative = Int64[]
    sp_iter = sp_term_iterator(rank, Val(N))
    highest_terms = Dict{semi_elementary_monomial{N},Int128}()

    for sp ∈ sp_iter
        empty!(highest_terms)
        for term ∈ p.s_terms[sp]
            highest_terms[term] = p.coeffs[term]
        end
        num_highest = N - findfirst(i -> i == sp[end], sp) + 1
        factor = ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], N)

        count_occurrences!(factor, f_keys, f_values)
        reset!(iter, f_values, num_highest)

        for way ∈ iter
            canonical_placement!(f_values, way, representative)
            new_term = ntuple(i -> factor[i] + representative[i], N)
            coeff = 1
            for i = 2:length(f_keys)
                if f_keys[i-1] == f_keys[i] - 1
                    coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                end
            end
            for (sp, c) in highest_terms
                push!(p, semi_elementary_monomial(new_term, sp.powers), -coeff * c)
            end
        end

        for (sp, c) in highest_terms
            tmp = ntuple(i -> i == num_highest ? sp.powers[i] + 1 : sp.powers[i], Val(N))
            push!(p, semi_elementary_monomial(factor, tmp), c)
        end

        pop!(p.s_terms, sp)
    end
    return p
end


function push_timeit!(p::semi_elementary_polynomial4{N}, x::semi_elementary_monomial{N}, v::Union{Integer,Rational}, timer) where {N}
    if !haskey(p.coeffs, x)
        @timeit timer "set coeff" p.coeffs[x] = v
        @timeit timer "alloc arr" (!haskey(p.s_terms, x.sp_term) && (p.s_terms[x.sp_term] = []))
        @timeit timer "push term" (x.sp_term[N] != 0 && push!(p.s_terms[x.sp_term], x))
    elseif p.coeffs[x] != -v
        @timeit timer "add coeff" p.coeffs[x] += v
    else
        @timeit timer "pop coeff" pop!(p.coeffs, x)
    end
end

function decompose13_timeit(x::semi_elementary_monomial{N}) where {N}
    rank = sum(x.sp_term)
    p = semi_elementary_polynomial4(N)
    x.sp_term == (0, 0, 0, 0, 0, 0, 40) && sizehint!(p.coeffs, 60000)
    push!(p, x, 1)
    timer = TimerOutput()

    @timeit timer "setup" begin
        iter = ways_place_containers_iterator([1], 1)
        f_keys = Int64[]
        f_values = Int64[]
        representative = Int64[]
        sp_iter = sp_term_iterator(rank, Val(N))
        highest_terms = Dict{semi_elementary_monomial{N},Int128}()
    end

    for sp ∈ sp_iter
        @timeit timer "get highest coeff" begin
            empty!(highest_terms)
            for term ∈ p.s_terms[sp]
                highest_terms[term] = p.coeffs[term]
            end
        end

        num_highest = N - findfirst(i -> i == sp[end], sp) + 1
        factor = @timeit timer "create factor" ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], Val(N))
        #@code_warntype ntuple(i -> i > N - num_highest ? sp[i] - 1 : sp[i], Val(N))
        #tmp = [factor...]

        @timeit timer "reset iter and count occurences" begin
            count_occurrences!(factor, f_keys, f_values)
            reset!(iter, f_values, num_highest)
        end
        println("f_keys:",f_keys," f_values:",f_values)

        for way ∈ iter
            @timeit timer "calc canonical placement" canonical_placement!(f_values, way, representative)
            #@code_warntype ntuple(i -> tmp[i] + representative[i], Val(N))
            #return
            new_term = @timeit timer "create new term" ntuple(i -> factor[i] + representative[i], Val(N))
            coeff = 1
            for i = 2:length(f_keys)
                if f_keys[i-1] == f_keys[i] - 1
                    coeff *= binomial(f_values[i] - way[i] + way[i-1], way[i-1])
                end
            end
            @timeit timer "push" begin
                for (sp, c) in highest_terms
                    push_timeit!(p, semi_elementary_monomial(new_term, sp.powers), -coeff * c, timer)
                end
            end
        end

        @timeit timer "push" begin
            for (sp, c) in highest_terms
                tmp = @timeit timer "create new term 2" ntuple(i -> i == num_highest ? sp.powers[i] + 1 : sp.powers[i], Val(N))
                push_timeit!(p, semi_elementary_monomial(factor, tmp), c, timer)
            end
        end

        @timeit timer "dealloc" pop!(p.s_terms, sp)
    end
    show(timer)
    return p
end

