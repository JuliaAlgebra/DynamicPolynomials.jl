using Test
import MultivariatePolynomials as MP

@testset "Variable and Monomial tests" begin
    @testset "Variable macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1] u[1:n, 1:n-1]
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
        @test zero_term(x) isa MP.Term{Int,monomial_type(XT)}
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
        @test_throws ErrorException DynamicPolynomials.Monomial(2)
        @test (@inferred DynamicPolynomials.Monomial(1)) isa monomial_type(XT)
        @test DynamicPolynomials.Monomial(1) == 1
        @test_throws ErrorException convert(monomial_type(XT), 2)
        @test (@inferred convert(monomial_type(XT), 1)) isa monomial_type(XT)
        @test convert(monomial_type(XT), 1) == 1
        @test_throws ErrorException convert(monomial_type(YT), 2)
        @test (@inferred convert(monomial_type(YT), 1)) isa monomial_type(YT)
        @test convert(monomial_type(YT), 1) == 1
        @polyvar x
        @test_throws ArgumentError monomial_type(XT)([x], [1, 0])
        @test zero_term(x^2) isa MP.Term{Int,monomial_type(x)}
        @test zero(x^2) isa polynomial_type(XT)
        @test one(x^2) isa monomial_type(XT)

        @polyvar y
        @test DynamicPolynomials.Monomial([x, y], [1, 0]) == x
        @test x != DynamicPolynomials.Monomial([x, y], [0, 1])
    end
    @testset "MonomialVector" begin
        @polyvar x y
        @test_throws AssertionError monomial_vector_type(x)([x], [[1], [1, 0]])
        X = MonomialVector([x, 1, x * y])
        @test variables(X) == [x, y]
        @test X.Z == [[0, 0], [1, 0], [1, 1]]
        @test monomial_vector_type(x)([1]) isa monomial_vector_type(x)
        @ncpolyvar nc
        @test monomial_vector_type(nc)([1]) isa monomial_vector_type(nc)
        a = [1, 5, 3]
        b, Y = monomial_vector(a, X)
        @test b === a
        @test Y === X
        σ, Y = sort_monomial_vector(X)
        @test σ == 1:length(X)
        @test Y === X
        XX = @inferred merge_monomial_vectors([X, X])
        @test XX isa monomial_vector_type(x)
        @test X == XX

        @test 2 != MonomialVector([x, y], 1)
        @test x != MonomialVector([x, y], 1)
        @test MonomialVector([x, y], [[0, 0], [1, 0]]) ==
              MonomialVector([x], [[0], [1]])
    end
    @testset "Non-commutative" begin
        @ncpolyvar x
        @test_throws ArgumentError monomial_type(x)([x], [1, 0])
        @test_throws AssertionError monomial_vector_type(x)([x], [[1], [1, 0]])
    end
    @testset "NC Variable * Monomial" begin
        @ncpolyvar x y z
        m = y * DynamicPolynomials.Monomial([y, z, x, z], [0, 0, 2, 1])
        @test variables(m) == [y, z, x, z]
        @test m.z == [1, 0, 2, 1]
        m = x * DynamicPolynomials.Monomial([z, y, y, z], [0, 0, 2, 1])
        @test variables(m) == [z, x, y, y, z]
        @test m.z == [0, 1, 0, 2, 1]
        m = x * DynamicPolynomials.Monomial([y, z, y, z], [0, 0, 2, 1])
        @test variables(m) == [y, z, x, y, z]
        @test m.z == [0, 0, 1, 2, 1]
    end
    @testset "NC Monomial * Variable" begin
        @ncpolyvar x y z
        m = DynamicPolynomials.Monomial([x, z, x, y], [2, 1, 0, 0]) * y
        @test variables(m) == [x, z, x, y]
        @test m.z == [2, 1, 0, 1]
        m = DynamicPolynomials.Monomial([x, y, y, x], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, y, z, x]
        @test m.z == [2, 1, 0, 1, 0]
        m = DynamicPolynomials.Monomial([x, y, x, y], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, z, x, y]
        @test m.z == [2, 1, 1, 0, 0]
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

        m = DynamicPolynomials.Monomial([x, y, z], [1, 2, 3])
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
        @test x == DynamicPolynomials.MP.map_exponents!(div, x^1, x * y^2)
        for z in [x * y, x^2, y^2]
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
    @testset "Noncommutative div" begin
        @ncpolyvar x y
	err = ErrorException("Not implemented yet")
	@test_throws err div(x, x * y)
    end
end
