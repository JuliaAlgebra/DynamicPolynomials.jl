@testset "Variable order" begin
    @polyvar x
    z = x
    @polyvar x
    @test z != x
end
@testset "Same variables grlex" begin
    @test DynamicPolynomials._exponents_compare([1, 0, 1], [1, 1, 0], MP.Graded{MP.LexOrder}) == -1
    @test DynamicPolynomials._exponents_compare([1, 1, 0], [1, 0, 1], MP.Graded{MP.LexOrder}) == 1
end
