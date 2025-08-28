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

@testset "Fast path cases: $ord" for ord in [
        InverseLexOrder,
        LexOrder,
        Graded{InverseLexOrder},
        Graded{LexOrder},
        Reverse{InverseLexOrder},
        Reverse{LexOrder},
        Graded{Reverse{InverseLexOrder}},
        Graded{Reverse{LexOrder}},
        Reverse{Graded{InverseLexOrder}},
        Reverse{Graded{LexOrder}},
    ]
    @polyvar x y z monomial_order=ord

    # allocation tests vary between Julia versions, so they're upper bounds
    @testset "Polynomial + constant" begin
        poly = 2 * x^2 + 3 * x * y + z * y^2
        poly2 = copy(poly)
        result = poly + 2
        MA.operate!(+, poly, 2)
        @test isequal(poly, result)
        # down from 576 using the generic method
        @test (@allocated MA.operate!(+, poly2, 2)) <= 272
        # subsequent additions don't allocate
        @test (@allocated MA.operate!(+, poly2, 2)) == 0

        # also test `-`
        poly = 2 * x^2 + 3 * x * y + z * y^2
        poly2 = copy(poly)
        result = poly - 2
        MA.operate!(-, poly, 2)
        @test isequal(poly, result)
        # down from 576 using the generic method
        @test (@allocated MA.operate!(-, poly2, 2)) <= 272
        # subsequent additions don't allocate
        @test (@allocated MA.operate!(-, poly2, 2)) == 0
    end

    @testset "Polynomial + Variable" begin
        poly = 2 * x^2 + 3 * x * y + z * y^2
        poly2 = copy(poly)
        result = poly + x
        MA.operate!(+, poly, x)
        @test isequal(poly, result)
        # down from 18752 using the generic method
        # 368 or 304 depending on ordering, more for different version
        # pre is especially bad
        @test (@allocated MA.operate!(+, poly2, x)) <= 400
        # down from 1904 using the generic method
        @test (@allocated MA.operate!(+, poly2, x)) <= 144

        # also test `-`
        poly = 2 * x^2 + 3 * x * y + z * y^2
        poly2 = copy(poly)
        result = poly - x
        MA.operate!(-, poly, x)
        @test isequal(poly, result)
        # down from 18752 using the generic method
        # 368 or 304 depending on ordering
        @test (@allocated MA.operate!(-, poly2, x)) <= 400
        # down from 1904 using the generic method
        @test (@allocated MA.operate!(-, poly2, x)) <= 144
    end
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
