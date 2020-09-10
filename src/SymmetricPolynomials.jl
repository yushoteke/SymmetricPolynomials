module SymmetricPolynomials

using Combinatorics
using TupleTools
using DataStructures

export symmetric_polynomial,decompose

include("elementary symmetric polynomial.jl")
include("elementary monomial.jl")
include("elementary polynomial.jl")
include("symmetric polynomial.jl")
include("monomial.jl")
include("utilities.jl")
include("decompose.jl")


end
