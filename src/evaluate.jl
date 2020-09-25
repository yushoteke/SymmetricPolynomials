
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

function evaluate_monomial(sp_values,m::semi_elementary_monomial)
    if !is_elementary(m)
        throw(ArgumentError("currently only support evaluating elementary monomials"))
    elseif length(sp_values) != dim(m)
        throw(DimensionMismatch("first argument must have same length"))
    end
    return prod([sp_values[i]^m.powers[i] for i=1:length(sp_values)])
end

function evaluate_polynomial(sp_values,m::semi_elementary_polynomial)
    return sum([evaluate_monomial(sp_values,k)*v for (k,v) in m.terms])
end
