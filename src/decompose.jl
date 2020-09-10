include("elementary polynomial.jl")
include("monomial.jl")
include("utilities.jl")

function simplify_sp_coeff_array(x)
    length(x) <= 1 && return x
    sorted = sort(x)
    tmp = eltype(x)[sorted[1]]
    for (i,x) in enumerate(sorted[2:end])
        if mergable(tmp[end][1],x[1])
            tmp[end] = (tmp[end][1],tmp[end][2]+x[2])
        else
            if tmp[end][2] == 0
                tmp[end] = x
            else
                push!(tmp,x)
            end
        end
    end
    return tmp
end

function multiply_two_elementary(x,y)
    x1,y1 = num_highest(x),num_highest(y)
    x1 < y1 && return multiply_two_elementary(y,x)
    N = dim(x)
    x0,y0 = N-x1,N-y1

    terms = Tuple{symmetric_polynomial{N},Rational{Int}}[]
    for num2 = max(y1-x0,0):y1
        num0 = x0 - (y1 - num2)
        tmp = ntuple(i->(if i>N-num2
                            return 2
                        elseif i>num0
                            return 1
                        else
                            return 0
                        end),N)
        #tmp = tuple(zeros(Int,num0)...,ones(Int,num1)...,2*ones(Int,num2)...)
        push!(terms,(symmetric_polynomial(tmp),binomial(y1+x1-2num2,y1-num2)))
    end
    return simplify_sp_coeff_array(terms)
end

#=
function multiply_one_elementary(x,y)
    #suppose y is elementary
    N = dim(x)
    terms = Tuple{symmetric_polynomial{N},Rational{Int64}}[]
    y1 = num_highest(y)
    y0 = N - y1
    c = symmetry_factor(x.exponents)
    orders = ordering_m1_n0(y1,y0)
    for i = 1:size(orders,1)
        sp = symmetric_polynomial(add_tuple_array(x.exponents,orders[i,:])...)
        coeff = symmetry_factor(sp.exponents)//c
        push!(terms,(sp,coeff))
    end
    return simplify_sp_coeff_array(terms)
end
=#

function multiply_one_elementary(x,y)
    N = dim(x)
    terms = Tuple{symmetric_polynomial{N},Rational{Int64}}[]
    y1 = num_highest(y)
    y0 = N - y1
    c = symmetry_factor(x.exponents)

    container_sizes = [1]
    for i=2:N
        if x.exponents[i]==x.exponents[i-1]
            container_sizes[end] += 1
        else
            push!(container_sizes,1)
        end
    end
    distribution_ways = ways_place_containers(container_sizes,y1)
    for way in distribution_ways
        representative = canonical_placement(container_sizes,way)
        sp = symmetric_polynomial((x.exponents.+representative)...)
        coeff = symmetry_factor(sp.exponents) * prod([binomial(container_sizes[i],way[i]) for i=1:length(way)]) // c
        push!(terms,(sp,coeff))
    end
    return simplify_sp_coeff_array(terms)
end

tables = Dict{Int64,Dict}()
for i=1:20
    tables[i] = Dict{symmetric_polynomial{i},elementary_polynomial{i}}()
end

decompose(x::symmetric_polynomial{N},solution_table=tables) where {N} = decompose_(x,tables[N])


function decompose_(x::symmetric_polynomial{N},solution_table:: Dict{symmetric_polynomial{N},elementary_polynomial{N}}) where {N}
    if is_elementary(x)
        tmp = elementary_symmetric_polynomial(N,sum(x.exponents))
        return convert(elementary_polynomial{N},tmp)
    end
    haskey(solution_table,x) && return solution_table[x]

    num_high = num_highest(x)
    factor1 = ntuple(i->i > N - num_high ? 1 : 0,N)
    factor2 = x.exponents .- factor1
    sp1 = elementary_symmetric_polynomial(N,num_high)
    sp2 = symmetric_polynomial(factor2...)
    initial_terms = is_elementary(sp2) ? multiply_two_elementary(sp2,sp1) :  multiply_one_elementary(sp2,sp1)
    solution = decompose_(sp2,solution_table) * sp1
    coef = 0
    for sp in initial_terms
        if !mergable(x,sp[1])
            solution -= sp[2] * decompose_(sp[1],solution_table)
        else
            coef = sp[2]
        end
    end
    solution /= coef
    solution_table[x] = solution
    return solution::elementary_polynomial{N}
end
