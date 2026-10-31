using Test
using DynamicPolynomials
import MultivariatePolynomials as MP
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
    @test (@inferred ordering(x[1])) === order()
    @test (@inferred ordering(typeof(x[1]))) === order()
    @test issorted(monomials(x[1], 0:2))
end

function _test_less(a, b)
    @test a < b
    @test b > a
end

function _test_monomials(vars, degs, exp)
    monos = @inferred monomials(vars, degs)
    @test monos == exp
    @test issorted(monos)
end

@testset "Noncommutative enumeration" begin
    @ncpolyvar x y
    _test_monomials((x, y), 0:2, [x^0, y, x, y * x, y^2, x * y, x^2])
    @ncpolyvar x y monomial_order = LexOrder
    @test_throws ArgumentError monomials((x, y), 2)
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

@testset "Noncommutative zero powers (DynamicPolynomials#118)" begin
    @ncpolyvar a x
    for (v, w) in ((a, x), (x, a))
        for p in (
            w^0 * v,
            v * w^0,
            (w * v)^0 * v,
            v * (w * v)^0,
            monomial([w, v, w], [0, 1, 0]),
        )
            @test (@inferred p == v)
            @test (@inferred v == p)
            @test (@inferred isequal(p, v))
            @test (@inferred isequal(v, p))
            @test p != w
            @test w != p
            @test hash(p) == hash(v)
            @test length(Set{Any}((p, v))) == 1
            @test Dict{Any,Int}(p => 1)[v] == 1
            @test Dict{Any,Int}(v => 1)[p] == 1
        end
    end
end

function _padded_noncommutative_word(m, vars)
    padded_vars = eltype(vars)[]
    padded_exps = Int[]
    for (v, e) in powers(m)
        push!(padded_vars, vars[end])
        push!(padded_exps, 0)
        if e >= 2
            append!(padded_vars, (v, vars[1], v))
            append!(padded_exps, (1, 0, e - 1))
        else
            push!(padded_vars, v)
            push!(padded_exps, e)
        end
    end
    push!(padded_vars, vars[1])
    push!(padded_exps, 0)
    return monomial(padded_vars, padded_exps)
end

@testset "Noncommutative word identity (DynamicPolynomials#131)" begin
    @ncpolyvar x[1:3]
    enumerated = monomials(x, 0:4)
    words = [[v for (v, e) in powers(m) for _ in 1:e] for m in enumerated]
    direct = [prod(word; init = one(x[1])) for word in words]
    padded = [_padded_noncommutative_word(m, x) for m in direct]

    @test issorted(enumerated)
    for i in eachindex(direct)
        a = direct[i]
        for b in (enumerated[i], padded[i])
            @test (@inferred a == b)
            @test (@inferred b == a)
            @test (@inferred isequal(a, b))
            @test iszero(@inferred cmp(a, b))
            @test iszero(@inferred cmp(b, a))
            for seed in (UInt(0), UInt(1), typemax(UInt))
                @test (@inferred hash(a, seed)) == hash(b, seed)
            end
            @test length(Set([a, b])) == 1
            @test Dict(a => 42)[b] == 42
        end
    end

    # The enumerated monomials share a padded variable sequence. Comparing
    # their exponent vectors gives an independent reference for the order.
    for i in eachindex(direct), j in eachindex(direct)
        a, b = direct[i], enumerated[j]
        expected = words[i] == words[j]
        @test (a == b) == expected
        @test isequal(a, b) == expected
        @test iszero(cmp(a, b)) == expected
        order = sign(cmp(ordering(a), exponents(enumerated[i]), exponents(b)))
        @test sign(cmp(a, b)) == order
        @test sign(cmp(padded[i], b)) == order
    end

    a = x[2] * x[1]
    b = monomial([x[1], x[2], x[3], x[1], x[2]], [0, 1, 0, 1, 0])
    p = 2a + 3b
    @test nterms(p) == 1
    @test coefficient(p, a) == coefficient(p, b) == 5
    @test iszero(p - 5a)

    a = x[1]^1000000 * x[2]
    b = monomial([x[3], x[1], x[2]], [0, 1000000, 1])
    @test (@inferred a == b)
    @test iszero(@inferred cmp(a, b))
    @test (@inferred hash(a)) == hash(b)
    @test collect(DynamicPolynomials._NoncommutativePowers(b)) ==
          [(x[1], 1000000), (x[2], 1)]
    @test isempty(collect(DynamicPolynomials._NoncommutativePowers(one(x[1]))))
    for m in (monomial([x[2]], [1]), monomial(x[1]))
        @test hash(m) == hash(variable(m))
    end
    @test hash(monomial([x[1], x[2], x[1]], [0, 1, 0])) == hash(x[2])
    @test hash(monomial([x[1], x[2]], [0, 0])) == hash(1)

    @ncpolyvar u v monomial_order = LexOrder
    a = u * v
    b = monomial([v, u, v, u], [0, 1, 1, 0])
    @test a == b
    @test iszero(cmp(a, b))
    @test cmp(u^1, v^2) > 0

    @ncpolyvar r s monomial_order = InverseLexOrder
    @test cmp(r * s, r^2) == cmp(InverseLexOrder(), [1, 1], [2, 0])
    @test cmp(r^2, r * s) < 0
end
