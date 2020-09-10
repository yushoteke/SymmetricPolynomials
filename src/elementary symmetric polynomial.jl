struct elementary_symmetric_polynomial{N}
    order::Int64
end

dim(x::elementary_symmetric_polynomial{N}) where {N} = N
elementary_symmetric_polynomial(N,M) = 0<=M<=N ? elementary_symmetric_polynomial{N}(M) : throw(DomainError("0 <= arg2 <= arg1 is violated"))
num_highest(x::elementary_symmetric_polynomial{N}) where {N} = x.order
