
"""
    ways_place_containers(x,m)

Given containers of sizes x, and m balls, enumerate all the ways to place the balls
into the containers without violating the constraints.

Lets illustrate this function through an example. Suppose we want to find
all the ways to place 4 balls into 4 containers with capacities [2,1,3,4]

We define f(m,e) as the number of solutions of placing m balls into the
e leftmost containers. Then we could build a table, vertically m and horizontally e.

The base cases are
m  e|   1   2   3   4
--------------------------
0   |   1   1   1   1
1   |   1
2   |   1
3   |   0
4   |   0

And then we could fill in the rest of the table using the recursive formula
f(m,e)=∑f(m-i,e-1)  where i=0:min(m,c[e]), where c[e] is the capacity of the e^th container.

m  e|   1   2   3   4
--------------------------
0   |   1   1   1   1
1   |   1   2   3   4
2   |   1   2   5   9
3   |   0   1   6   15
4   |   0   0   5   20

Then we could use this table to build the actual solutions with another builder function.
"""
function ways_place_containers(x, m)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with
    #positive capacity x_i. Now we want to place m balls into these containers
    #without breaking the capacity
    #this solution will use a dynamic programming approach.
    #First, find the number of solutions, then allocate an array, then
    #build that table using the solution table and a builder function.

    #base cases
    solution_table = zeros(Int64, m + 1, length(x))
    solution_table[1, 1:end] .= 1
    for r = 2:m+1
        solution_table[r, 1] = r - 1 <= x[1]
    end

    #build solution table
    for c = 2:length(x)
        for r = 2:m+1
            solution_table[r, c] = sum(solution_table[r-min(r - 1, x[c]):r, c-1])
        end
    end

    ways = zeros(Int64, solution_table[end, end], length(x))
    table_builder!(ways, solution_table, x, 1, length(x), m)
    return ways
end

function table_builder!(output, solution_table, container_sizes, row, col, balls_left)
    (balls_left <= 0) && return
    r1 = row
    for balls_to_place = 0:min(balls_left, container_sizes[col])
        if col == 1
            output[row, col] = balls_left
            continue
        end
        ways = solution_table[balls_left+1-balls_to_place, col-1]
        ways <= 0 && continue
        output[r1:r1+ways-1, col] .= balls_to_place
        table_builder!(output, solution_table, container_sizes, r1, col - 1, balls_left - balls_to_place)

        r1 = r1 + ways
    end
end

mutable struct ways_place_containers_iterator
    capacities::Vector{Int64}
    num_balls::Int64
    vec::Vector{Int64}
end

function ways_place_containers_iterator(x, m)
    tmp = zeros(eltype(x), length(x))
    return ways_place_containers_iterator(x, m, tmp)
end

"""
    Here is the iteration strategy, suppose x=[2,1,3,4,5],m=11
    suppose we have some configuration, say 21341, then we can imagine this as a number, 
    2x^4 + x^3 + 3x^2 + 4x + 1. For x big enough, this creates an ordering, where the larger it is on the left side, the larger it is

    Then, we first find the largest configuration, which is 21341 in our case, and we successively find the "next smallest number" 
    by lowering the right most lowerable column. 

    column 4 is the largest lowerable column
    21341
    21332
    21323
    21314
    21305

    Now, column 3 is the largest lowerable column, we lower column 3 by 1, and fill in the rest greedily
    21242

    And now column 4 is the largest lowerable column again, and we repeat this process
    21233
    21224
    21215
    21143
    21134
    21125
    21044
    ........

    A column is lowerable if its value is not 0 and the larger columns are not saturated.
    When there are no lowerable columns, the iterator ends

    This implementation simply returns the vector that the iterator is holding. I understand that this is bad practice,
    because it is extremely unsafe. However, I haven't found a "return const reference" type of solution, and this is one 
    way I find that doesn't allocate any new memory
"""
function Base.iterate(x::ways_place_containers_iterator, state=0)

    # if state = 0, then that means the array needs to be initialized
    if state == 0
        sum(x.capacities) < x.num_balls && return nothing
        fill_in_greedily!(x.vec, x.capacities, 1, x.num_balls)
        return (x.vec, state + 1)
    end

    #if vec has already been initialized, then find the rightmost column that is lowerable
    cum_capacity = x.capacities[end]
    cum_balls = x.vec[end]
    for i = length(x.vec)-1:-1:1
        if x.vec[i] != 0 && cum_balls != cum_capacity
            x.vec[i] -= 1
            fill_in_greedily!(x.vec, x.capacities, i + 1, cum_balls + 1)
            return (x.vec, state + 1)
        else
            cum_capacity += x.capacities[i]
            cum_balls += x.vec[i]
        end
    end

    #if reach this point, then that means there are no more lowerable column, and the iterator ends
    return nothing
end

function fill_in_greedily!(tgt_vec, constraints, start_ind, num_balls)
    for i = start_ind:length(tgt_vec)
        tgt_vec[i] = min(constraints[i], num_balls)
        num_balls -= tgt_vec[i]
    end
end

function reset!(x::ways_place_containers_iterator, constraints, balls)
    empty!(x.capacities)
    for c ∈ constraints
        push!(x.capacities, c)
    end
    x.num_balls = balls
    empty!(x.vec)
    resize!(x.vec, length(constraints))
end

"""
    canonical_placement(x,y)

Given containers of size x, and number of balls to place into each container y,
return the placement as a tuple of 1s and 0s.

Example:

canonical_placement([2,3,4],[1,1,1]) returns

(1,0,1,0,0,1,0,0,0)

Which one could think of as 10|100|1000

"""
function canonical_placement(x, y)
    #assume x=[x_1,x_2,x_3...x_n] respresents n containers each with positive capacity x_i.
    #assume y=[y_1,y_2,y_3...y_n] means placed 0<=y_i<=x_i balls into container i
    #returns one way this could have happened.
    #For example, x=[1,2,6],y=[1,1,2] could be 1|10|010100 or 1|01|110000
    tmp = zeros(Int, sum(x))
    ind = 1
    for i = 1:length(x)
        tmp[ind:ind+y[i]-1] .= 1
        ind += x[i]
    end
    return tmp
end

function canonical_placement!(x, y, tgt_vec)
    empty!(tgt_vec)
    for (xx, yy) ∈ zip(x, y)
        for _ = 1:yy
            push!(tgt_vec, 1)
        end
        for _ = 1:xx-yy
            push!(tgt_vec, 0)
        end
    end
end

function count_occurrences(factor)
    #assume factor is a sorted tuple like (0,0,0,1,1,3,3)
    #returns two arrays, the first one keys, the second one
    #number of occurrences, like for this case,
    #[0,1,3],[3,2,2]
    x = [factor[1]]
    y = [1]
    for i = 2:length(factor)
        if factor[i] == factor[i-1]
            y[end] += 1
        else
            push!(x, factor[i])
            push!(y, 1)
        end
    end
    return x, y
end

function count_occurrences!(factor, x, y)
    empty!(x)
    empty!(y)
    push!(x, factor[1])
    push!(y, 1)
    for i = 2:length(factor)
        if factor[i] == factor[i-1]
            y[end] += 1
        else
            push!(x, factor[i])
            push!(y, 1)
        end
    end
end