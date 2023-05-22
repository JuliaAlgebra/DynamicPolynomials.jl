@testset "Variable order" begin
    @polyvar x
    z = x
    @polyvar x
    @test z != x
end
@testset "README example" begin
    p(x, y, z) = sprint(show, MIME"text/plain"(), 4x*y^2*z + 4z^2 - 5x^3 + 7x^2*z^2)
    @polyvar x y z monomial_order = LexOrder
    @test p(x, y, z) == "4z² + 4xy²z + 7x²z² - 5x³"
    @polyvar x y z
    @test p(x, y, z) == "4z² - 5x³ + 4xy²z + 7x²z²"
    @polyvar x y z monomial_order = Graded{Reverse{InverseLexOrder}}
    @test p(x, y, z) == "4z² - 5x³ + 7x²z² + 4xy²z"
end
