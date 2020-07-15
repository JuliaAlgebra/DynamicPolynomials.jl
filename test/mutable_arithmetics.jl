using Test

import MutableArithmetics
const MA = MutableArithmetics

using DynamicPolynomials

@testset "Mutable Arithmetics with $T" for T in [Int, BigInt]
    @polyvar x y z
    p = T(2)x^2 + T(4)y^2
    q = T(3)x * y
    @test p === MA.add!(p, q)
    @test q == T(3)x * y
    @test nterms(p) == 3
    @test issorted(monomials(p), rev=true)
    @test p == T(2)x^2 + T(3)x * y + T(4)y^2
    @test p === MA.add!(p, x * y)
    @test nterms(p) == 3
    @test issorted(monomials(p), rev=true)
    @test p == T(2)x^2 + T(4)x * y + T(4)y^2
    @test p === MA.add!(p, T(2)x^3 + T(3)x^2 + T(5)x + T(2))
    @test nterms(p) == 6
    @test issorted(monomials(p), rev=true)
    @test p == T(2)x^3 + T(5)x^2 + T(4)x * y + T(4)y^2 + T(5)x + T(2)
    q = T(5)z^3 + T(1)x*z + T(4)z^2 + T(3)z
    @test p === MA.add!(p, q)
    @test nterms(p) == 10
    @test issorted(monomials(p), rev=true)
    @test p == T(2)x^3 + T(5)x^2 + T(4)x * y + T(4)y^2 + T(5)x + T(2) + q

    @testset "Issue #62" begin
        @polyvar x y
        p = x^2 + x + 1
        q = rem(p, [x^2-y])
        @test q == x + y + 1 
    end
end


