
"""
    semi_elementary_monomial(x,y)

First, given ``n`` variables, we say ``e_i`` is the elementary symmetric polynomial or order i with coefficient 1.\n
For example, when ``n=3``,

``e_1 = x+y+z``

``e_2 = xy+xz+yz``

``e_3 = xyz``

When ``n=4``,

``e_1 = w+x+y+z``

``e_2 = wx+wy+wz+xy+xz+yz``

``e_3 = wxy+wxz+wyz+xyz``

``e_4 = wxyz``

and so on

From the users perspective, usually the ``n`` is obvious, so it is omitted.

Then, ``S(a_1,a_2...a_n)`` is the symmetric polynomial in ``n`` variables with powers ``0 \\leq a_1 \\leq a_2...\\leq a_n`` and coefficient 1\n
For example,

``S(0,1,2) = x^2y + x^2z + y^2z``

``S(0,0,1,2) = w^2x + w^2y + w^2z + x^2y + x^2z + y^2z``

``S(0,0,0,0,5) = v^5 + w^5 + x^5 + y^5 + z^5``

and so on

A semi elementary monomial is the product of\n
    1. A non elementary symmetric polynomial with power 1
    2. An elementary monomial

For example, the follow are semi elementary

``S(0,0,3)e_1^3e_2^1e_3^0``

``S(0,0,1,2)e_1^3e_2^0e_3^1e_4^2``

These would be represented as two tuples,

``S(0,0,3)e_1^3e_2^1e_3^0           => (0,0,3),(3,1,0)``

``S(0,0,1,2)e_1^3e_2^0e_3^1e_4^2  => (0,0,1,2),(3,0,1,2)``

This is where the name semi comes from, in that nearly all terms in the monomial are elementary except 1.

When the user wants to create a semi elementary monomial, he or she needs to provide two tuples with non-negative integers.
The first tuple will automatically be sorted by the constructor, the second will be provided as it is.

    semi_elementary_monomial((0,2,1),(3,4,5))

returns the same object as

    semi_elementary_monomial((0,1,2),(3,4,5))

Besides sorting, some trivial simplifications will also be made.
First, ``e_n`` terms will be automatically pulled out, for example

    semi_elementary_monomial((2,3,4),(0,0,0))

will return

    semi_elementary_monomial((0,1,2),(0,0,2))

Which is basically saying

``(x^2y^3z^4 + x^2z^3y^4 + y^2x^3z^4 + y^2z^3x^4 + z^2x^3y^4 + z^2x^4y^3) = (x^2y + x^2z + y^2z) * (xyz)^2``

Second, if the first term is already elementary, it will be merged into the second term, for example

    semi_elementary_monomial((0,1,1),(0,0,0))

will return

    semi_elementary_monomial((0,0,0),(0,1,0))
"""
struct semi_elementary_monomial{N}
    sp_term::NTuple{N,Int64}
    powers::NTuple{N,Int64}

    function semi_elementary_monomial(x::NTuple{N,Int64},y::NTuple{N,Int64}) where {N}
        t = TupleTools.sort(x)
        if t[end] - t[1] == 1
            c = sum(t) - N * t[1]
            tmp2 = ntuple(N) do i
                    if i==c
                        return y[i] + 1
                    elseif i==N
                        return y[i] + t[1]
                    else
                        return y[i]
                    end
                end
            tmp1 = ntuple(i->0,N)
        elseif t[1] > 0
            tmp2 = ntuple(i->i==N ? y[i]+t[1] : y[i],N)
            tmp1 = t .- t[1]
        else
            tmp1,tmp2 = t,y
        end
        return new{N}(tmp1,tmp2)
    end
end


dim(x::semi_elementary_monomial{N}) where {N} = N
is_elementary(x::semi_elementary_monomial) = x.sp_term[end]==0
Base.:(==)(x::semi_elementary_monomial,y::semi_elementary_monomial) = (x.sp_term==y.sp_term) && (x.powers==y.powers)
function Base.isless(x::semi_elementary_monomial{N},y::semi_elementary_monomial{N}) where {N}
    for i=N:-1:1
        x.sp_term[i] != y.sp_term[i] && return x.sp_term[i] < y.sp_term[i]
    end
    s1,s2 = sum(x.powers),sum(y.powers)
    return s1!=s2 ? s1 < s2 : x.powers < y.powers
    #return x.powers < y.powers
end

function to_string(x::semi_elementary_monomial)
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
Base.show(io::IO, ::MIME{Symbol("text/plain")}, x::semi_elementary_monomial) = join(io,to_string(x))
