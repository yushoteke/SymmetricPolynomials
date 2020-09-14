using SymmetricPolynomials
using Test
import SymmetricPolynomials: semi_elementary_polynomial,push!


@testset "Decompose Elementary Symmetric Polynomials" begin
    include("decomposeTests.jl")
end
