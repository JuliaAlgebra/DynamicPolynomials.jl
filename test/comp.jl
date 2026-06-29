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
function _less(a, b)
    @test a < b
    @test b > a
    @test cmp(monomial(a), b) < 0
    @test cmp(b, monomial(a)) > 0
end
@testset "Issue 152" begin
    @polyvar x y monomial_order=LexOrder
    _less(x, x * y)
    _less(x * y^2, x^2)
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

function _test_less(a, b)
    @test a < b
    @test b > a
end

function _test_monomials(vars, degs, exp)
    # Without `collect`, `exp` is promoted to a `MonomialVector`
    # which sorts it so it doesn't test the order
    @test collect(monomials(vars, degs)) == exp
end

@testset "LexOrder" begin
    @polyvar x y monomial_order = LexOrder
    _test_less(y, y^2)
    _test_less(x^0, y)
    _test_less(y^2, x)
    _test_less(x * y^2, x^2)
    _test_monomials([x, y], 3, [y^3, y^2 * x, y * x^2, x^3])
    _test_monomials([x, y], 1:2, [y, y^2, x, x * y, x^2])
    _test_monomials([x, y], [0, 1, 3], [1, y, y^3, x, y^2 * x, y * x^2, x^3])
end

@testset "InverseLexOrder" begin
    @polyvar x y monomial_order = InverseLexOrder
    _test_less(y, y^2)
    _test_less(x^0, y)
    _test_less(x, y^2)
    _test_less(x^2, x * y^2)
    _test_monomials([x, y], 3, [x^3, x^2 * y, x * y^2, y^3])
end

@testset "Graded{Reverse{LexOrder}}" begin
    @polyvar x y monomial_order = Graded{Reverse{LexOrder}}
    _test_less(y, y^2)
    _test_less(x^0, y)
    _test_less(x, y^2)
    _test_less(x^2, x * y^2)
    _test_monomials([x, y], 3, [x^3, x^2 * y, x * y^2, y^3])
    _test_monomials([x, y], 1:2, [x, y, x^2, x * y, y^2])
    _test_monomials([x, y], [0, 1, 3], [1, x, y, x^3, x^2 * y, x * y^2, y^3])
end

@testset "Graded{Reverse{InverseLexOrder}}" begin
    @polyvar x y monomial_order = Graded{Reverse{InverseLexOrder}}
    _test_less(y, y^2)
    _test_less(x^0, y)
    _test_less(x, y^2)
    _test_less(x^2, x * y^2)
    _test_monomials([x, y], 3, [y^3, y^2 * x, y * x^2, x^3])
    _test_monomials([x, y], 1:2, [y, x, y^2, x * y, x^2])
    _test_monomials([x, y], [0, 1, 3], [1, y, x, y^3, x*y^2, x^2*y, x^3])
end

@testset "Comparison with NaN" begin
    @polyvar p
    poly1 = NaN + p
    poly2 = NaN + p
    @test poly1 != poly2
    @test (@allocated poly1 != poly2) == 0
    @test poly1 != poly1
    @test (@allocated poly1 != poly1) == 0
    @test isequal(poly1, poly2)
    @test (@allocated isequal(poly1, poly2)) == 0
    @test isequal(poly1, poly1)
    @test (@allocated isequal(poly1, poly1)) == 0
end
