
function ordering_m1_n0(m,n)
    arr = zeros(Int,binomial(m+n,m),m+n)
    set_1s_0s(arr,1,1,m+n,m,n)
    return arr
end
function set_1s_0s(arr,h,w,wmax,m,n)
    x = binomial(m+n-1,n)
    arr[h:h+x-1,w] .= 1
    if w!=wmax
        set_1s_0s(arr,h,w+1,wmax,m-1,n)
        set_1s_0s(arr,h+x,w+1,wmax,m,n-1)
    end
end

function ways_place_containers(x,m)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with
    #positive capacity x_i. Now we want to place m balls into these containers
    #without breaking the capacity

    #base case, assume x <= sum(x), so there is at least 1 valid solution
    length(x)==1 && return Array{Int,1}[[m]]
    results = Array{Int,1}[]
    maximum = min(x[end],m)
    minimum = max(0,m - sum(x[1:end-1]))
    for balls_in_last_container = minimum:maximum
        subsolutions = ways_place_containers(x[1:end-1],m-balls_in_last_container)
        if subsolutions!= nothing
            for i in subsolutions
                push!(i,balls_in_last_container)
            end
        end
        append!(results,subsolutions)
    end
    return results
end

function canonical_placement(x,y)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with positive capacity x_i.
    #assume y=[y_1,y_2,y_3...y_n] means placed 0<=y_i<=x_i balls into container i
    #returns one way this could have happened.
    #For example, x=[1,2,6],y=[1,1,2] could be 1|10|010100 or 1|01|110000
    tmp = zeros(Int,sum(x))
    starting_locations = cumsum([1,x...])
    for (i,x) in enumerate(starting_locations[1:end-1])
        tmp[x:x-1+y[i]] .= 1
    end
    return tmp
end

add_tuple_array(x::NTuple{N,T},y) where {N,T<:Any} = ntuple(i->x[i]+y[i],N)
sorttuple(x::NTuple) = length(x)<10 ? TupleTools.sort(x) : tuple(sort([x...])...)

function get_valid_exponents(dim,rank)
    tmp = integer_partitions(rank)
    too_long = findall(x->length(x)>dim,tmp)
    deleteat!(tmp,too_long)
    for i in tmp
        if length(i) < dim
            append!(i,zeros(Int64,dim-length(i)))
        end
        reverse!(i)
    end
    return tmp
end

function count_occurrences(x)
    #count number of each item in a sorted array
    #for example, [3,3,3,3,3,4,4,5] -> [5,2,1]
    container_sizes = [1]
    for i=2:length(x)
        if x[i]==x[i-1]
            container_sizes[end] += 1
        else
            push!(container_sizes,1)
        end
    end
    return container_sizes
end

function symmetry_factor(exponents)
    occur = count_occurrences(exponents)
    return prod(factorial.(occur))
end

#converts an term to elementary representation
#for example, (0,1,1) -> (0,1,0),(0,0,1) -> (1,0,0), (1,1,1)->(0,0,1)
(elementary_representation(x::NTuple{N,T})::NTuple{N,T}) where {N,T}  = ntuple(i->i==sum(x) ? one(T) : zero(T),N)
