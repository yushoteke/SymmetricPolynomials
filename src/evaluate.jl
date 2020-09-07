#elementary symmetric polynomial values
esp_values = Dict{Tuple{Int,Int},Float64}()
#first entry is dimension of elementary symmetric polynomial
#second entry is order of rank of elementary symmetric polynomial

function evaluate(x::symmetric_polynomial)
    if is_elementary(x)
        k = (dim(x),sum(x.exponents))
        if haskey(esp_values,k)
            return esp_values[k]
        else
            throw(ErrorException("user hasn't defined ",k))
        end
    else
        p = decompose_symmetric_polynomial(x)
        return evaluate(p)
    end
end

evaluate(x::polynomial) = sum(evaluate.(x.terms))
evaluate(x::monomial) = x.coeff * prod([evaluate(i[1])^i[2] for i in x.factors])

function sp_value_from_primitive_values(x...)
    #=
    assume the values provided are x1,x2...
    generates the list of values for e1,e2...
    For example, x=1,y=2,z=3
    Then e1=x+y+z = 6,e2=xy+yz+xz = 11, e3=xyz = 6
    =#
    length(x) < 1 && return []
    p = [1,x[1]]
    for i = 2:length(x)
        tmp = [0,x[i]*p...]
        tmp[1:length(p)] += p
        p=tmp
    end
    return p[2:end]
end

function add_primitive_values(x...)
    N = length(x)
    v = sp_value_from_primitive_values(x...)
    for i=1:N
        esp_values[(N,i)] = v[i]
    end
end
