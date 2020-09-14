
@testset "Basic polynomial functionalities" begin
    @testset "convert monomial to polynomial" begin
        x = semi_elementary_monomial((1,2,3),(0,0,0))
        p = to_polynomial(x)
        s = semi_elementary_polynomial(3)
        push!(s,semi_elementary_monomial((0,1,2),(0,0,1)),1)
        @test s == p
    end
end
