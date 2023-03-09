using DataStructures
"""
  A structure which specializes in dictionaries of the form
  Dict{X,Vector{Y}}

  Internally, keeps a pool of arrays. When a key is popped,
  the corresponding value(which is a vector) is returned to the 
  pool. When a new key is created, it is assigned a vector from 
  the pool. 

  Internally, the arrays are also assigned tiers. Each tier's array 
  capacity is 4 times the previous tier's array capacity. When an array 
  of tier X(A_X) needs to resize, first checks if there's an empty tier X+1 array 
  (A_{X+1}) in the pool. If there is, then contents of A_X are copied into 
  A_{X+1}, A_X is emptied and returned to the pool. Or else, A_X itself is
  resized and promoted into tier X+1. 


"""
struct ArrayPoolDict

end

"""
  the above is the ultimate goal, but for now, just use a simpler strategy
"""
struct ArrayPoolDict1{X,Y}
  dict::Dict{X,Vector{Y}}
  pool::Stack{Vector{Y}}
end

ArrayPoolDict1{X,Y}() where {X,Y} = ArrayPoolDict1{X,Y}(Dict(), Stack{Vector{Y}}())

function Base.push!(apd::ArrayPoolDict1{X,Y}, key::X, value::Y) where {X,Y}
  if !haskey(apd.dict, key)
    arr = borrow!(apd)
    apd.dict[key] = arr
  end
  push!(apd.dict[key], value)
end

function Base.getindex(apd::ArrayPoolDict1{X,Y}, key::X) where {X,Y}
  return apd.dict[key]
end

function borrow!(apd::ArrayPoolDict1{X,Y}) where {X,Y}
  if isempty(apd.pool)
    return Y[]
  else
    return pop!(apd.pool)
  end
end

function Base.pop!(apd::ArrayPoolDict1{X,Y}, key::X) where {X,Y}
  arr = apd.dict[key]
  empty!(arr)
  push!(apd.pool, arr)
  delete!(apd.dict, key)
end

