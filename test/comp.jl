@testset "Variable order" begin
    @polyvar x
    z = x
    @polyvar x
    @test z != x
end
@testset "Same variables grlex" begin
    @test DynamicPolynomials._samevars_grlex([1, 0, 1], [1, 1, 0]) == -1
    @test DynamicPolynomials._samevars_grlex([1, 1, 0], [1, 0, 1]) == 1
end
