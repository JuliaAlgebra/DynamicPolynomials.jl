using Test
using DynamicPolynomials
import DynamicPolynomials: Commutative, CreationOrder # to test hygiene
@testset "Variable order" begin
    @polyvar x
    z = x
    @polyvar x
    @test z != x
    order = Commutative{CreationOrder}
    @polyvar x variable_order = order
    @test z != x
end
@testset "README example" begin
    function p(x, y, z)
        return sprint(
            show,
            MIME"text/plain"(),
            4x * y^2 * z + 4z^2 - 5x^3 + 7x^2 * z^2,
        )
    end
    @polyvar x y z monomial_order = LexOrder
    @test p(x, y, z) == "4z² + 4xy²z + 7x²z² - 5x³"
    @polyvar x y z
    @test p(x, y, z) == "4z² - 5x³ + 4xy²z + 7x²z²"
    @polyvar x y z monomial_order = Graded{Reverse{InverseLexOrder}}
    @test p(x, y, z) == "4z² - 5x³ + 7x²z² + 4xy²z"
end
# See https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/138
# Also tests `ordering`
@testset "InverseLexOrder" begin
    order = Graded{InverseLexOrder}
    @polyvar x[1:2] monomial_order = order
    @test ordering(x[1]) == order
    @test issorted(monomials(x[1], 0:2))
end
