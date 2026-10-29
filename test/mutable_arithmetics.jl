using Test

import MutableArithmetics
const MA = MutableArithmetics
import MultivariatePolynomials as MP

using DynamicPolynomials

@testset "Mutable Arithmetics with $T" for T in [Int, BigInt]
    @polyvar x y z
    p = T(2)x^2 + T(4)y^2
    q = T(3)x * y
    @test p === MA.add!!(p, q)
    @test q == T(3)x * y
    @test nterms(p) == 3
    @test issorted(monomials(p))
    @test p == T(2)x^2 + T(3)x * y + T(4)y^2
    @test p === MA.add!!(p, x * y)
    @test nterms(p) == 3
    @test issorted(monomials(p))
    @test p == T(2)x^2 + T(4)x * y + T(4)y^2
    @test p === MA.add!!(p, T(2)x^3 + T(3)x^2 + T(5)x + T(2))
    @test nterms(p) == 6
    @test issorted(monomials(p))
    @test p == T(2)x^3 + T(5)x^2 + T(4)x * y + T(4)y^2 + T(5)x + T(2)
    q = T(5)z^3 + T(1)x * z + T(4)z^2 + T(3)z
    @test p === MA.add!!(p, q)
    @test nterms(p) == 10
    @test issorted(monomials(p))
    @test p == T(2)x^3 + T(5)x^2 + T(4)x * y + T(4)y^2 + T(5)x + T(2) + q

    @testset "Issue #62" begin
        @polyvar x y
        p = x^2 + x + 1
        q = rem(p, [x^2 - y])
        @test q == x + y + 1
    end
end

@testset "In-place addition and subtraction: $ord" for ord in [
    InverseLexOrder,
    LexOrder,
    Graded{InverseLexOrder},
    Graded{LexOrder},
    Graded{Reverse{InverseLexOrder}},
    Graded{Reverse{LexOrder}},
]
    @polyvar x y z monomial_order=ord

    # Arithmetic now uses StarAlgebras instead of the former DP fast paths.
    # Check insertion and subsequent updates without assuming their allocations.
    @testset "Polynomial $op $name" for op in (+, -),
        (name, operand) in (("constant", 2), ("variable", x))

        poly = 2 * x^2 + 3 * x * y + z * y^2
        result = op(poly, operand)
        @test MA.operate!(op, poly, operand) === poly
        @test isequal(poly, result)

        result = op(result, operand)
        @test MA.operate!(op, poly, operand) === poly
        @test isequal(poly, result)
    end
end

@testset "Invalid monomial ordering: $ord" for ord in [
    Reverse{InverseLexOrder},
    Reverse{LexOrder},
    Reverse{Graded{InverseLexOrder}},
    Reverse{Graded{LexOrder}},
]
    @polyvar x monomial_order = ord
    @test_throws ArgumentError 2x
end

@testset "Non-concrete in-place polynomial addition" begin
    @polyvar p q r s
    p1 = MP.polynomial(r - q, Number) + MP.polynomial(3//5 * p^2, Number)
    p2 = MP.polynomial(1//2 + s, Number) + MP.polynomial(p^2, Number)
    result = p1 + p2
    @test isequal(result, MA.operate!(+, p1, p2))

    p1 = MP.polynomial(r - q, Number) + MP.polynomial(3//5 * p^2, Number)
    p2 = MP.polynomial(1//2 + s, Number) + MP.polynomial(p^2, Number)
    @test isequal(result, MA.operate!(+, p2, p1))
end
