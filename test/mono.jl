using Test
using DynamicPolynomials
import MultivariatePolynomials as MP

@testset "Variable and Monomial tests" begin
    @testset "Variable macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:(n-1)] u[1:n, 1:(n-1)]
        VT = typeof(y)
        @test x isa Vector{VT}
        @test y isa VT
        @test z isa Vector{VT}
        @test u isa Matrix{VT}
        @test length(x) == 3
        @test length(z) == 2
        @test size(u) == (3, 2)
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
        @test u[1, 1] > u[2, 1] > u[2, 2]
        dummy = VT("dummy")
        @test dummy isa VT
        @test name(dummy) == "dummy"

        @polyvar a[1:5, 1:3, 1:2]
        @test size(a) == (5, 3, 2)
    end
    @testset "Variable macro tuple return" begin
        vars = @polyvar x y z
        @test vars isa Tuple
        @test vars == (x, y, z)

        vars = @ncpolyvar x y z
        @test vars isa Tuple
        @test vars == (x, y, z)
    end

    @testset "variable_union_type" begin
        @polyvar x
        XT = typeof(x)
        @test DynamicPolynomials.MP.variable_union_type(x) == XT
        @test DynamicPolynomials.MP.variable_union_type(x^2) == XT
        @test DynamicPolynomials.MP.variable_union_type(2x) == XT
        @test DynamicPolynomials.MP.variable_union_type(x + 1) == XT
        @ncpolyvar y
        YT = typeof(y)
        @test DynamicPolynomials.MP.variable_union_type(y) == YT
        @test DynamicPolynomials.MP.variable_union_type(y^2) == YT
        @test DynamicPolynomials.MP.variable_union_type(2y) == YT
        @test DynamicPolynomials.MP.variable_union_type(y + 1) == YT
    end
    @testset "Variable" begin
        @polyvar x
        XT = typeof(x)
        @test zero_term(XT) == 0
        @test zero(XT) == 0
        @ncpolyvar y
        YT = typeof(y)
        @test one(YT) == 1
        @polyvar x
        @test zero_term(x) isa term_type(XT)
        @test zero(x) isa polynomial_type(XT)
        @test one(x) isa monomial_type(XT)
    end
    @testset "Monomial" begin
        @ncpolyvar ncp
        YT = typeof(ncp)
        @test zero_term(YT) == 0
        @test zero(YT) == 0
        @polyvar cp
        XT = typeof(cp)
        @test one(XT) == 1
        for VT in (XT, YT)
            MT = monomial_type(VT)
            @test (@inferred constant_monomial(VT)) isa MT
            @test one(MT) == 1
            @test_throws InexactError convert(MT, 2)
            @test (@inferred convert(MT, 1)) isa MT
            @test convert(MT, 1) == 1
        end
        @polyvar x
        @test_throws ArgumentError monomial([x], [1, 0])
        @test zero_term(x^2) isa term_type(x)
        @test zero(x^2) isa polynomial_type(XT)
        @test one(x^2) isa monomial_type(XT)

        @polyvar y
        @test monomial([x, y], [1, 0]) == x
        @test x != monomial([x, y], [0, 1])
    end
    @testset "MonomialVector" begin
        @polyvar x y
        @test_throws ArgumentError monomial_vector([
            monomial([x], e) for e in ([1], [1, 0])
        ])
        X = monomial_vector([x, 1, x * y])
        @test variables(X) == [x, y]
        @test [exponents(m, [x, y]) for m in X] == [[0, 0], [1, 0], [1, 1]]
        @test monomial_vector_type(x)([1]) isa monomial_vector_type(x)
        @ncpolyvar nc
        @test monomial_vector_type(nc)([1]) isa monomial_vector_type(nc)
        a = [1, 5, 3]
        b, Y = monomial_vector(a, X)
        @test b == a
        @test Y == X
        σ, Y = sort_monomial_vector(X)
        @test σ == 1:length(X)
        @test Y == X
        XX = @inferred merge_monomial_vectors([X, X])
        @test XX isa monomial_vector_type(x)
        @test X == XX

        @test 2 != monomials([x, y], 1)
        @test x != monomials([x, y], 1)
        @test monomial_vector([
            monomial([x, y], e) for e in ([0, 0], [1, 0])
        ]) == monomial_vector([monomial([x], e) for e in ([0], [1])])
    end
    @testset "Non-commutative" begin
        @ncpolyvar x
        @test_throws ArgumentError monomial([x], [1, 0])
        @test_throws ArgumentError monomial_vector([
            monomial([x], e) for e in ([1], [1, 0])
        ])
    end
    @testset "NC Variable * Monomial" begin
        @ncpolyvar x y z
        m = y * monomial([y, z, x, z], [0, 0, 2, 1])
        @test variables(m) == [y, x, z]
        @test exponents(m) == [1, 2, 1]
        m = x * monomial([z, y, y, z], [0, 0, 2, 1])
        @test variables(m) == [x, y, z]
        @test exponents(m) == [1, 2, 1]
        m = x * monomial([y, z, y, z], [0, 0, 2, 1])
        @test variables(m) == [x, y, z]
        @test exponents(m) == [1, 2, 1]
    end
    @testset "NC Monomial * Variable" begin
        @ncpolyvar x y z
        m = monomial([x, z, x, y], [2, 1, 0, 0]) * y
        @test variables(m) == [x, z, y]
        @test exponents(m) == [2, 1, 1]
        m = monomial([x, y, y, x], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, z]
        @test exponents(m) == [2, 1, 1]
        m = MP.monomial([x, y, x, y], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, z]
        @test exponents(m) == [2, 1, 1]
    end

    @testset "Antidifferentiation" begin
        @ncpolyvar x y z

        m = x
        mi = DynamicPolynomials.MP.antidifferentiate(m, y)
        @test mi == x * y

        # Antidifferentiation is product => Integral coefficients
        @test MP.coefficient_type(mi) == Int

        # General antidifferentiation => Rational coefficients
        m = x^3
        mi = DynamicPolynomials.MP.antidifferentiate(m, x)
        @test mi == (x^4 / 4)
        @test MP.coefficient_type(mi) == Rational{Int}

        m = monomial([x, y, z], [1, 2, 3])
        mi = DynamicPolynomials.MP.antidifferentiate(m, z)
        @test mi == (x*y^2*z^4) / 4
        @test MP.coefficient_type(mi) == Rational{Int}
    end

    @testset "Evaluation" begin
        @polyvar x y
        @test (x^2 * y)(3, 2) == 18
        @test (x^2 * y)((3, 2)) == 18
        @test (x^2 * y)([3, 2]) == 18
        @test (x^2)(3) == 9
        @test (x)(3) == 3
    end
    @testset "TODO remove when added to MP" begin
        @polyvar x y
        @test x == DynamicPolynomials.MP.map_exponents!(div, x * y^0, x * y^2)
        @test_throws ArgumentError MP.map_exponents!(div, x^1, x * y^2)
        for z in [x * y, x^2 * y^0, x^0 * y^2]
            @test y ==
                  DynamicPolynomials.MP.map_exponents_to!(z, -, x * y^2, x * y)
        end
    end
    # TODO add to MP
    @testset "Indexing with boolean" begin
        @polyvar x y
        X = monomials([x, y], 2)
        @test X[[true, false, true]] == monomial_vector([x^2, y^2])
        X = monomials([x, y], 0:1)
        @test filter(mono -> degree(mono) == 1, X) == monomial_vector([x, y])
        @test filter(mono -> degree(mono) == 0, X) == monomial_vector([x^0])
    end
    @testset "NC exponents with variables" begin
        @ncpolyvar x y z
        # Monomial with subset of variables, exponents mapped to full variable list
        m = monomial([x, z], [2, 3])
        e = MP.exponents(m, [x, y, z])
        @test e == [2, 0, 3]
        @test length(e) == 3

        # Single variable monomial mapped to larger variable set
        m = monomial([y], [4])
        e = MP.exponents(m, [x, y, z])
        @test e == [0, 4, 0]

        # Monomial using all variables
        m = monomial([x, y, z], [1, 2, 3])
        e = MP.exponents(m, [x, y, z])
        @test e == [1, 2, 3]

        # Constant monomial (no variables)
        m = monomial(typeof(x)[], Int[])
        e = MP.exponents(m, [x, y, z])
        @test e == [0, 0, 0]

        # Noncommutative monomial with repeated variables
        m = monomial([x, y, x], [1, 2, 3])
        e = MP.exponents(m, [x, y, x])
        @test e == [1, 2, 3]
    end
    @testset "Noncommutative div" begin
        @ncpolyvar x y
        err = ErrorException("Not implemented yet")
        @test_throws err div(x, x * y)
    end
end
