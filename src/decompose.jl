
function multiply_two_elementary(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N}
    x1,y1 = sum(x.exponents),sum(y.exponents)
    x1 < y1 && return multiply_two_elementary(y,x)
    x0,y0 = N-x1,N-y1

    terms = monomial{N}[]
    for num2 = max(y1-x0,0):y1
        c = binomial(y1+x1-2num2,y1-num2)
        num0 = x0 - (y1 - num2)
        tmp = zeros(Int64,dim(x))
        tmp[num0+1:end] .= 1
        tmp[end-num2+1:end] .= 2
        push!(terms,monomial(c,symmetric_polynomial(tmp)))
    end
    return polynomial(terms)
end

function multiply_one_elementary(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N}
    #suppose y is elementary
    terms = monomial{N}[]
    y1 = sum(y.exponents)
    y0 = N - y1
    c = symmetry_factor(x.exponents)
    expo = [x.exponents...]
    for order in ordering_m1_n0(y1,y0)
        tmp = symmetric_polynomial(expo+order)
        push!(terms,to_monomial(tmp) * (symmetry_factor(tmp.exponents)//c))
    end
    return polynomial(simplify_monomial_array(terms))
end

decompose_table = Dict{symmetric_polynomial,polynomial}()

decompose(x::polynomial) =sum([decompose(i) for i in x.terms])
decompose(x::monomial) = x.coeff * prod([decompose(i[1])^i[2] for i in x.factors])

function decompose(x::symmetric_polynomial)
    is_elementary(x) && return to_polynomial(x)
    haskey(decompose_table,x) && return decompose_table[x]

    expo = [x.exponents...]
    num_same_as_largest = length(expo) - findfirst(x->x==expo[end],expo) + 1
    factor1 = zeros(Int64,length(expo))
    factor1[end-num_same_as_largest+1:end] .= 1
    factor2 = expo - factor1

    sp1 = symmetric_polynomial(factor1)
    sp2 = symmetric_polynomial(factor2)
    initial_guess = is_elementary(sp2) ? multiply_two_elementary(sp2,sp1) :  multiply_one_elementary(sp2,sp1)

    G = to_polynomial(to_monomial(sp1)*to_monomial(sp2))

    ind = findfirst(t->summable(to_monomial(x),t),initial_guess.terms)
    similar_term = initial_guess.terms[ind]
    c = similar_term.coeff

    rewrite_x = (G-initial_guess+to_polynomial(similar_term))*(1/c)
    solution = decompose(rewrite_x)
    decompose_table[x] = solution
    return solution
end
