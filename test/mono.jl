@testset "PolyVar and Monomial tests" begin
    @testset "polyvar macro index set" begin
        n = 3
        @polyvar x[1:n] y z[1:n-1]
        @test x isa Vector{PolyVar{true}}
        @test y isa PolyVar{true}
        @test z isa Vector{PolyVar{true}}
        @test length(x) == 3
        @test length(z) == 2
        @test x[1] > x[2] > x[3] > y > z[1] > z[2]
    end
    @testset "PolyVar" begin
        @test zeroterm(PolyVar{true}) == 0
        @test zero(PolyVar{true}) == 0
        @test one(PolyVar{false}) == 1
        @polyvar x
        @test zeroterm(x) isa Term{true, Int}
        @test zero(x) isa Polynomial{true, Int}
        @test one(x) isa Monomial{true}
    end
    @testset "Monomial" begin
        @test zeroterm(Monomial{false}) == 0
        @test zero(Monomial{false}) == 0
        @test one(Monomial{true}) == 1
        @polyvar x
        @test_throws ArgumentError Monomial{true}([x], [1,0])
        @test zeroterm(x^2) isa Term{true, Int}
        @test zero(x^2) isa Polynomial{true, Int}
        @test one(x^2) isa Monomial{true}
    end
    @testset "MonomialVector" begin
        @polyvar x y
        @test_throws ArgumentError MonomialVector{true}([x], [[1], [1,0]])
        X = MonomialVector([x, 1, x*y])
        @test variables(X) == [x, y]
        @test X.Z == [[1, 1], [1, 0], [0, 0]]
        @test MonomialVector{true}([1]) isa MonomialVector{true}
        @test MonomialVector{false}([1]) isa MonomialVector{false}
        a = [1, 5, 3]
        b, Y = monovec(a, X)
        @test b === a
        @test Y === X
        σ, Y = sortmonovec(X)
        @test σ == 1:length(X)
        @test Y === X
        @test mergemonovec([X, X]) == X
    end
    @testset "Non-commutative" begin
        @ncpolyvar x
        @test_throws ArgumentError Monomial{false}([x], [1,0])
        @test_throws ArgumentError MonomialVector{false}([x], [[1], [1,0]])
    end
    @testset "NC PolyVar * Monomial" begin
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
    @testset "NC Monomial * PolyVar" begin
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
end
