
"""
A semi elementary monomial is the product of
1.a non elementary symmetric polynomial with power 1
2.an elementary monomial

For example 3 * S(2,1,0,0) * e1^3 * e2 * e3 * e4^2

The non elementary term is always sorted
"""

struct semi_elementary_monomial{N}
    sp_term::NTuple{N,Int64}
    powers::NTuple{N,Int64}

    function semi_elementary_monomial(x,y)
        x = sorttuple(x)
        N = length(x)
        if x[1]!=0
            #automatically factors out e_n terms like xyz relative to x+y+z and xy+xz+yz
            #after that, check if the remaining is elementary
            x,y = x .- x[1],ntuple(i->i==N ? x[1]+y[i] : y[i],N)
            if x[end]==1
                x,y =ntuple(i->0,N), y .+ elementary_representation(x)
            end
        end
        return new{N}(x,y)
    end
end

dim(x::semi_elementary_monomial{N}) where {N} = N

is_elementary(x::semi_elementary_monomial) = x.sp_term[end]==0

#assumes y is elementary, returns the product
mul_unsafe(x::semi_elementary_monomial{N},y::semi_elementary_monomial{N}) where {N} = semi_elementary_monomial(x.sp_term,x.powers.+y.powers)

function Base.isless(x::semi_elementary_monomial{N},y::semi_elementary_monomial{N}) where {N}
    for i=N:-1:1
        x.sp_term[i] != y.sp_term[i] && return x.sp_term[i] < y.sp_term[i]
    end
    return x.powers < y.powers
end

function semi_elementary_monomial_to_string(x::semi_elementary_monomial)
    if !is_elementary(x)
        head = "S"*string(x.sp_term)
    else
        head = ""
    end
    for (i,n) in enumerate(x.powers)
        if n == 0
            continue
        elseif n == 1
            head *= "e$i "
        else
            head *= "e$i^$n "
        end
    end
    return head
end
Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::semi_elementary_monomial) = join(io,semi_elementary_monomial_to_string(x))
