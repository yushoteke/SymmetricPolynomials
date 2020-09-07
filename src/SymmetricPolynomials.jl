module SymmetricPolynomials

using Combinatorics
using TupleTools
using DataStructures

export symmetric_polynomial,decompose,to_polynomial

include("symmetric polynomial.jl")
include("monomial.jl")
include("polynomial.jl")
include("utilities.jl")
include("decompose.jl")
include("evaluate.jl")


end
