using SymmetricPolynomials
using Test
import SymmetricPolynomials: semi_elementary_polynomial, dim, is_elementary, to_string
import SymmetricPolynomials: push!

@testset "Monomials" begin
    include("monomialTests.jl")
end

@testset "Decompose Elementary Symmetric Polynomials" begin
    include("decomposeTests.jl")
end
