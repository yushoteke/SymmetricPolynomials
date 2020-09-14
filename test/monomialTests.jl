
@testset "Basic Monomial Functionalities" begin
    @testset "invalid arguments 1" begin
        @test_throws DomainError semi_elementary_monomial((-1,0,0),(0,0,0))
        @test_throws DomainError semi_elementary_monomial((0,0,0),(-1,0,0))
        @test_throws DimensionMismatch semi_elementary_monomial((0,0,0),(0,0,0,0))
        @test_throws DimensionMismatch semi_elementary_monomial((0,0,0,0),(0,0,0))
    end

    @testset "automatic sorting" begin
        x = semi_elementary_monomial((3,2,1),(1,0,0))
        @test x.sp_term == (0,1,2)
        x = semi_elementary_monomial((2,1,0),(1,0,0))
        @test x.sp_term == (0,1,2)
    end

    @testset "pull out common factors" begin
        x = semi_elementary_monomial((3,4,5),(0,0,0))
        @test x == semi_elementary_monomial((0,1,2),(0,0,3))
    end

    @testset "test <, 0 terms same" begin
        x = semi_elementary_monomial((3,5,5),(0,0,0))
        y = semi_elementary_monomial((3,4,6),(0,0,0))
        @test (x < y) == true
    end

    @testset "test <, 1 term same" begin
        x = semi_elementary_monomial((3,3,6),(0,0,0))
        y = semi_elementary_monomial((3,4,6),(0,0,0))
        @test (x < y) == true
    end

    @testset "test <, all sp terms same" begin
        x = semi_elementary_monomial((3,3,6),(0,0,0))
        y = semi_elementary_monomial((3,3,6),(1,0,0))
        @test (x < y) == true
    end

    @testset "test <, autosimplifying" begin
        x = semi_elementary_monomial((2,2,2),(1,0,0))
        y = semi_elementary_monomial((3,3,3),(0,0,0))
        @test (x < y) == false
    end

    @testset "test is_elementary" begin
        x = semi_elementary_monomial((3,3,3),(0,0,0))
        @test is_elementary(x) == true
        x = semi_elementary_monomial((1,0,0),(0,0,0))
        @test is_elementary(x) ==false
        x = semi_elementary_monomial((0,0,1),(0,0,0))
        @test is_elementary(x) ==false
    end

    @testset "test dim" begin
        x = semi_elementary_monomial((3,3,3),(0,0,0))
        @test dim(x) == 3
        x = semi_elementary_monomial((3,3,3,3),(0,0,0,0))
        @test dim(x) == 4
    end

    @testset "string conversion" begin
        x = semi_elementary_monomial((3,3,3),(0,0,0))
        @test to_string(x) == "e3^3 "
        x = semi_elementary_monomial((3,3,0),(0,0,0))
        @test to_string(x) == "S(0, 3, 3)"
        x = semi_elementary_monomial((3,3,1),(1,0,0))
        @test to_string(x) == "S(0, 2, 2)e1 e3 "
    end
end
