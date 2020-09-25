using SymmetricPolynomials
using Test
import SymmetricPolynomials: semi_elementary_polynomial, dim, is_elementary, to_string
import SymmetricPolynomials: push!, to_polynomial

@testset "Monomials" begin
    include("monomialTests.jl")
end

@testset "Polynomials" begin
    include("polynomialTests.jl")
end

@testset "Decompose" begin
    include("decomposeTests.jl")
end

@testset "Evaluate" begin
    include("evaluateTests.jl")
end
