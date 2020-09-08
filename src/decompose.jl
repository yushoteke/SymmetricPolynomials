
function multiply_two_elementary(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N}
    x1,y1 = sum(x.exponents),sum(y.exponents)
    x1 < y1 && return multiply_two_elementary(y,x)
    x0,y0 = N-x1,N-y1

    terms = monomial{N}[]
    for num2 = max(y1-x0,0):y1
        c = binomial(y1+x1-2num2,y1-num2)
        num0 = x0 - (y1 - num2)
        tmp = ntuple(i->(if i>N-num2
                            return 2
                        elseif i>num0
                            return 1
                        else
                            return 0
                        end),N)
        push!(terms,monomial(c,symmetric_polynomial(tmp...)))
    end
    return polynomial(terms)
end

function multiply_one_elementary(x::symmetric_polynomial{N},y::symmetric_polynomial{N}) where {N}
    #suppose y is elementary
    terms = monomial{N}[]
    y1 = sum(y.exponents)
    y0 = N - y1
    c = symmetry_factor(x.exponents)
    for order in ordering_m1_n0(y1,y0)
        tmp = symmetric_polynomial(add_tuple_array(x.exponents,order)...)
        #push!(terms,to_monomial(tmp) * (symmetry_factor(tmp.exponents)//c))
        push!(terms,monomial(symmetry_factor(tmp.exponents)//c,tmp))
    end
    return polynomial(simplify_monomial_array(terms))
end

tables = Dict{Int64,Dict}()
for i=1:20
    tables[i] = Dict{symmetric_polynomial{i},polynomial{i}}()
end

function decompose_symmetric_polynomial(x::symmetric_polynomial{N},decompose_table:: Dict{symmetric_polynomial{N},polynomial{N}}) where {N}
    is_elementary(x) && return to_polynomial(x)
    haskey(decompose_table,x) && return decompose_table[x]
    num_same_as_largest = N - findfirst(i->i==x.exponents[end],x.exponents) + 1
    factor1 = ntuple(i->i > N - num_same_as_largest ? 1 : 0,N)
    factor2 = x.exponents .- factor1
    sp1 = symmetric_polynomial(factor1...)
    sp2 = symmetric_polynomial(factor2...)
    initial_guess = is_elementary(sp2) ? multiply_two_elementary(sp2,sp1) :  multiply_one_elementary(sp2,sp1)
    G = to_polynomial(to_monomial(sp1)*to_monomial(sp2))
    ind = findfirst(t->summable(to_monomial(x),t),initial_guess.terms)
    similar_term = initial_guess.terms[ind]
    c = similar_term.coeff
    rewrite_x = (G-initial_guess+to_polynomial(similar_term))*(1/c)
    solution = decompose_polynomial(rewrite_x,decompose_table)
    decompose_table[x] = solution
    return solution
end

decompose_polynomial(x::polynomial{N},decompose_table:: Dict{symmetric_polynomial{N},polynomial{N}}) where {N} = sum(polynomial{N}[decompose_monomial(i,decompose_table) for i in x.terms])
decompose_monomial(x::monomial{N},decompose_table:: Dict{symmetric_polynomial{N},polynomial{N}}) where {N} = x.coeff * prod(polynomial{N}[decompose_symmetric_polynomial(i[1],decompose_table)^i[2] for i in x.factors])
