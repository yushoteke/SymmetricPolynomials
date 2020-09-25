
"""
    semi_elementary_monomial(x,y)

First, given n variables, we say e_i is the elementary symmetric polynomial or order i with coefficient 1.\n
For example,\n
    n=3, e1 = x+y+z,    e2 = xy+xz+yz,          e3 = xyz
    n=4, e1 = w+x+y+z,  e2 = wx+wy+wz+xy+xz+yz, e3 = wxy+wxz+wyz+xyz, e4 = wxyz
    ...
    and so on

Second, S(a1,a2...an) is the symmetric polynomial in n variables with powers 0 <= a1 <= a2...<= an and coefficient 1\n
For example,\n
    S(0,1,2) = x^2y + x^2z + y^2z
    S(0,0,1,2) = w^2x + w^2y + w^2z + x^2y + x^2z + y^2z
    S(0,0,0,0,4) = v^5 + w^5 + x^5 + y^5 + z^5
    ...
    and so on

A semi elementary monomial is the product of\n
    1.a non elementary symmetric polynomial with power 1
    2.an elementary monomial

For example,\n
    S(0,0,3) * e1^3 * e2^1 * e3^0            = (x^3 + y^3 + z^3) * (x + y + z)^3 * (xy + xz + yz)
    S(0,0,1,2) * e1^3 * e2^0 * e3^1 * e4^2   = (w^2x + w^2y + w^2z + x^2y + x^2z + y^2z) * (w + x + y + z)^3 * (wxy + wxz + wyz + xyz) * (wxyz)^2

These would be represented as two tuples, for example\n
    S(0,0,3) * e1^3 * e2^1 * e3^0           => (0,0,3),(3,1,0)
    S(0,0,1,2) * e1^3 * e2^0 * e3^1 * e4^2  => (0,0,1,2),(3,0,1,2)

This is where the name semi comes from, in that nearly all terms in the monomial are elementary except 1.

When the user wants to create a semi elementary monomial, he or she needs to provide two tuples with non-negative integers.
The first tuple will automatically be sorted by the constructor, the second will be provided as it is.

    semi_elementary_monomial((0,2,1),(3,4,5))

returns the same object as

    semi_elementary_monomial((0,1,2),(3,4,5))

Besides sorting, some trivial simplifications will also be made.
First, e_n terms will be automatically pulled out, for example

    semi_elementary_monomial((2,3,4),(0,0,0))

will return

    semi_elementary_monomial((0,1,2),(0,0,2))

Which is basically saying
    (x^2y^3z^4 + x^2z^3y^4 + y^2x^3z^4 + y^2z^3x^4 + z^2x^3y^4 + z^2x^4y^3) = (x^2y + x^2z + y^2z) * (xyz)^2

Second, if the first term is already elementary, it will be merged into the second term, for example

    semi_elementary_monomial((0,1,1),(0,0,0))

will return

    semi_elementary_monomial((0,0,0),(0,1,0))
"""
struct semi_elementary_monomial{N}
    sp_term::NTuple{N,Int64}
    powers::NTuple{N,Int64}

    function semi_elementary_monomial(x,y)
        (length(x) != length(y)) && throw(DimensionMismatch("Input tuples should have same length"))
        #(any(i->i<0,x) || any(i->i<0,y)) && throw(DomainError("all entries must be greater than 0"))
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
