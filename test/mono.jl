using Test

@testset "Variable and Monomial tests" begin
    @testset "Variable macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1] u[1:n, 1:n-1]
        @test x isa Vector{Variable{true}}
        @test y isa Variable{true}
        @test z isa Vector{Variable{true}}
        @test u isa Matrix{Variable{true}}
        @test length(x) == 3
        @test length(z) == 2
        @test size(u) == (3, 2)
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
        @test u[1, 1] > u[2, 1] > u[2, 2]

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
        @test DynamicPolynomials.MP.variable_union_type(x) == Variable{true}
        @test DynamicPolynomials.MP.variable_union_type(x^2) == Variable{true}
        @test DynamicPolynomials.MP.variable_union_type(2x) == Variable{true}
        @test DynamicPolynomials.MP.variable_union_type(x + 1) == Variable{true}
        @ncpolyvar y
        @test DynamicPolynomials.MP.variable_union_type(y) == Variable{false}
        @test DynamicPolynomials.MP.variable_union_type(y^2) == Variable{false}
        @test DynamicPolynomials.MP.variable_union_type(2y) == Variable{false}
        @test DynamicPolynomials.MP.variable_union_type(y + 1) ==
              Variable{false}
    end
    @testset "Variable" begin
        @test zero_term(Variable{true}) == 0
        @test zero(Variable{true}) == 0
        @test one(Variable{false}) == 1
        @polyvar x
        @test zero_term(x) isa Term{true,Int}
        @test zero(x) isa Polynomial{true,Int}
        @test one(x) isa Monomial{true}
    end
    @testset "Monomial" begin
        @test zero_term(Monomial{false}) == 0
        @test zero(Monomial{false}) == 0
        @test one(Monomial{true}) == 1
        if VERSION ≥ v"0.7-"
            @test_throws ErrorException Monomial(2)
            @test (@inferred Monomial(1)) isa Monomial{true}
            @test Monomial(1) == 1
            @test_throws ErrorException convert(Monomial{true}, 2)
            @test (@inferred convert(Monomial{true}, 1)) isa Monomial{true}
            @test convert(Monomial{true}, 1) == 1
            @test_throws ErrorException convert(Monomial{false}, 2)
            @test (@inferred convert(Monomial{false}, 1)) isa Monomial{false}
            @test convert(Monomial{false}, 1) == 1
        end
        @polyvar x
        @test_throws ArgumentError Monomial{true}([x], [1, 0])
        @test zero_term(x^2) isa Term{true,Int}
        @test zero(x^2) isa Polynomial{true,Int}
        @test one(x^2) isa Monomial{true}

        @polyvar y
        @test Monomial([x, y], [1, 0]) == x
        @test x != Monomial([x, y], [0, 1])
    end
    @testset "MonomialVector" begin
        @polyvar x y
        @test_throws AssertionError MonomialVector{true}([x], [[1], [1, 0]])
        X = MonomialVector([x, 1, x * y])
        @test variables(X) == [x, y]
        @test X.Z == [[0, 0], [1, 0], [1, 1]]
        @test MonomialVector{true}([1]) isa MonomialVector{true}
        @test MonomialVector{false}([1]) isa MonomialVector{false}
        a = [1, 5, 3]
        b, Y = monomial_vector(a, X)
        @test b === a
        @test Y === X
        σ, Y = sort_monomial_vector(X)
        @test σ == 1:length(X)
        @test Y === X
        XX = @inferred merge_monomial_vectors([X, X])
        @test XX isa MonomialVector{true}
        @test X == XX

        @test 2 != MonomialVector([x, y], 1)
        @test x != MonomialVector([x, y], 1)
        @test MonomialVector([x, y], [[0, 0], [1, 0]]) ==
              MonomialVector([x], [[0], [1]])
    end
    @testset "Non-commutative" begin
        @ncpolyvar x
        @test_throws ArgumentError Monomial{false}([x], [1, 0])
        @test_throws AssertionError MonomialVector{false}([x], [[1], [1, 0]])
    end
    @testset "NC Variable * Monomial" begin
        @ncpolyvar x y z
        m = y * Monomial([y, z, x, z], [0, 0, 2, 1])
        @test variables(m) == [y, z, x, z]
        @test m.z == [1, 0, 2, 1]
        m = x * Monomial([z, y, y, z], [0, 0, 2, 1])
        @test variables(m) == [z, x, y, y, z]
        @test m.z == [0, 1, 0, 2, 1]
        m = x * Monomial([y, z, y, z], [0, 0, 2, 1])
        @test variables(m) == [y, z, x, y, z]
        @test m.z == [0, 0, 1, 2, 1]
    end
    @testset "NC Monomial * Variable" begin
        @ncpolyvar x y z
        m = Monomial([x, z, x, y], [2, 1, 0, 0]) * y
        @test variables(m) == [x, z, x, y]
        @test m.z == [2, 1, 0, 1]
        m = Monomial([x, y, y, x], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, y, z, x]
        @test m.z == [2, 1, 0, 1, 0]
        m = Monomial([x, y, x, y], [2, 1, 0, 0]) * z
        @test variables(m) == [x, y, z, x, y]
        @test m.z == [2, 1, 1, 0, 0]
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
end
