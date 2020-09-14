
@testset "decompose" begin
    @testset "x^2+y^2+z^2" begin
        S = semi_elementary_monomial
        x = S((2,0,0),(0,0,0))
        poly = semi_elementary_polynomial(3)
        push!(poly,S((0,0,0),(2,0,0)),1)
        push!(poly,S((0,0,0),(0,1,0)),-2)
        @test decompose(x) == poly
    end

    @testset "x^3+y^3+z^3" begin
        S = semi_elementary_monomial
        x = S((3,0,0),(0,0,0))
        poly = semi_elementary_polynomial(3)
        push!(poly,S((0,0,0),(3,0,0)),1)
        push!(poly,S((0,0,0),(1,1,0)),-3)
        push!(poly,S((0,0,0),(0,0,1)),3)
        @test decompose(x) == poly
    end

    @testset "x^4+y^4+z^4" begin
        S = semi_elementary_monomial
        x = S((4,0,0),(0,0,0))
        poly = semi_elementary_polynomial(3)
        push!(poly,S((0,0,0),(4,0,0)),1)
        push!(poly,S((0,0,0),(2,1,0)),-4)
        push!(poly,S((0,0,0),(1,0,1)),4)
        push!(poly,S((0,0,0),(0,2,0)),2)
        @test decompose(x) == poly
    end

    @testset "x^5+y^5+z^5" begin
        S = semi_elementary_monomial
        x = S((5,0,0),(0,0,0))
        poly = semi_elementary_polynomial(3)
        push!(poly,S((0,0,0),(5,0,0)),1)
        push!(poly,S((0,0,0),(3,1,0)),-5)
        push!(poly,S((0,0,0),(2,0,1)),5)
        push!(poly,S((0,0,0),(1,2,0)),5)
        push!(poly,S((0,0,0),(0,1,1)),-5)
        @test decompose(x) == poly
    end

    @testset "w^5+x^5+y^5+z^5" begin
        S = semi_elementary_monomial
        x = S((5,0,0,0),(0,0,0,0))
        poly = semi_elementary_polynomial(4)
        push!(poly,S((0,0,0,0),(5,0,0,0)),1)
        push!(poly,S((0,0,0,0),(3,1,0,0)),-5)
        push!(poly,S((0,0,0,0),(2,0,1,0)),5)
        push!(poly,S((0,0,0,0),(1,2,0,0)),5)
        push!(poly,S((0,0,0,0),(1,0,0,1)),-5)
        push!(poly,S((0,0,0,0),(0,1,1,0)),-5)
        @test decompose(x) == poly
    end
end
