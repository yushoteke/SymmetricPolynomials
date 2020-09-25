
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
