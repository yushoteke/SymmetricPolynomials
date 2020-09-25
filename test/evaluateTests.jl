FP = SymmetricPolynomials.sp_value_from_primitive_values
EM = SymmetricPolynomials.evaluate_monomial
EP = SymmetricPolynomials.evaluate_polynomial

@testset "Evaluate" begin
    @testset "evaluate monomials" begin
        xyz = [1.0,2.0,3.0]
        sp_values = FP(xyz...)
        x = semi_elementary_monomial((0,0,0),(1,0,0))
        @test EM(sp_values,x) ≈ 6.0 atol=0.001
        x = semi_elementary_monomial((0,0,0),(0,1,0))
        @test EM(sp_values,x) ≈ 11.0 atol=0.001
        x = semi_elementary_monomial((0,0,0),(0,0,1))
        @test EM(sp_values,x) ≈ 6.0 atol=0.001
        x = semi_elementary_monomial((0,0,0),(1,2,3))
        @test EM(sp_values,x) ≈ 156816.0 atol=0.001
    end

    @testset "evaluate polynomials" begin
        xyz = [1.0,2.0,3.0]
        sp_values = FP(xyz...)
        x = semi_elementary_monomial((0,0,0),(1,0,0))
        y = semi_elementary_monomial((0,0,0),(0,1,0))
        p = to_polynomial(x)
        push!(p,y,1)
        @test EP(sp_values,p) ≈ 17.0 atol=0.001
        push!(p,y,-1)
        @test EP(sp_values,p) ≈ 6.0 atol=0.001
        z = semi_elementary_monomial((0,0,0),(1,2,3))
        push!(p,z,2)
        @test EP(sp_values,p) ≈ 2*156816+6 atol=0.001
    end

    @testset "evaluate small decomposed polynomials" begin
        xyz = [1.0,2.0,3.0]
        sp_values = FP(xyz...)
        x = semi_elementary_monomial((0,0,3),(0,0,0))
        p = decompose(x)
        @test EP(sp_values,p) ≈ 1^3+2^3+3^3 atol=0.001
        x = semi_elementary_monomial((0,0,5),(0,0,0))
        p = decompose(x)
        @test EP(sp_values,p) ≈ 1^5+2^5+3^5 atol=0.001
    end

    @testset "evaluate big decomposed polynomials" begin

        x = semi_elementary_monomial((0,0,0,0,0,0,0,0,0,25),(0,0,0,0,0,0,0,0,0,0))
        p = decompose(x)
        for i=1:30
            xyz = collect(range(BigInt(1)//1,stop=1+i//20,length=10))
            sp_values = FP(xyz...)
            expected = sum(xyz.^25)
            @test EP(sp_values,p) == expected
        end
    end
end
