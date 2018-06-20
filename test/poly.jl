@testset "Term and Polynomial tests" begin
    @testset "Term" begin
        @polyvar x
#       @test coefficienttype(1x) == Int
#       @test coefficienttype(1.0x^2) == Float64
#       @test coefficienttype(Term{true, Int}) == Int
        @test zeroterm(Term{false, Int}).α == 0
        @test one(Term{true, Int}).α == 1
        @polyvar x
        @test typeof(Term(1x)) == Term{true, Int}
        @test Term(1x) == 1x
        @test typeof(convert(Any, 1x)) == Term{true, Int}
        @test typeof(one(1x)) == Term{true, Int}
        @test typeof(zeroterm(1x)) == Term{true, Int}
        @test typeof(zero(1x)) == Polynomial{true, Int}
        @test typeof(one(1.0x)) == Term{true, Float64}
        @test typeof(zeroterm(1.0x)) == Term{true, Float64}
        @test typeof(zero(1.0x)) == Polynomial{true, Float64}

        @test typeof(constantterm(1, x)) == Term{true, Int}
        @inferred constantterm(1, x)
        @test typeof(polynomial(1:2, monomials([x], 1:2))) == Polynomial{true, Int}
    end

    @testset "Polynomial" begin
        #@test eltype(Polynomial{true, Int}) == Int
        @polyvar x
        @test_throws ArgumentError Polynomial{true, Int}([1, 2], [x])
        @test_throws ArgumentError Polynomial{true, Int}([1, 2], MonomialVector([x]))
        @test_throws InexactError Polynomial{true, Int}([1.5], [x])
        @test Polynomial(1 + x) == 1 + x
        @test typeof(one(1 + x)) == Polynomial{true, Int}
        @test typeof(zeroterm(1 + x)) == Term{true, Int}
        @test typeof(zero(1 + x)) == Polynomial{true, Int}
        @test typeof(one(1.0 + x)) == Polynomial{true, Float64}
        @test typeof(zeroterm(1.0 + x)) == Term{true, Float64}
        @test typeof(zero(1.0 + x)) == Polynomial{true, Float64}

        p = Polynomial([4, 9], [x, x*x])
        p.a == [9, 4]
        p.x[1] == x^2
        p.x[2] == x

        @inferred Polynomial(i -> float(i), [x, x*x])
        @inferred Polynomial(i -> 3 - float(i), MonomialVector([x*x, x]))
        for p in (Polynomial(i -> float(i), [x, x*x]),
                  Polynomial(i -> 3 - float(i), MonomialVector([x*x, x])))
            @test typeof(p) == Polynomial{true, Float64}
            @test p.a == [2.0, 1.0]
            @test p.x == MonomialVector([x^2, x])
        end

        @ncpolyvar ncpolyvar u v
        @inferred Polynomial(i -> i, [u, u*u, 1])
        p = Polynomial(i -> i, [u, u*u, 1])
        @test typeof(p) == Polynomial{false, Int}
        @test p.a == [2, 1, 3]
        @test p.x == MonomialVector([u^2, u, 1])

        @test u + v*u + 1 != v*u + u
        @test removemonomials(u + v*u + 1, [1, u*v]) == v*u + u
        @test removemonomials(u + u*v + 1, [u*v]) == 1 + u

        @inferred polynomial(2u)
        @inferred polynomial(2.0u, Int)
    end

    @testset "Evaluation" begin
        @polyvar x y
        @test (x^2+y^3)(2, 3) == 31
        @test (x^2+y^3)((2, 3)) == 31
        @test (x^2+y^3)([2, 3]) == 31
        @test (2x^2*y)(3,2) == 36
        @test (2x^2*y)((3,2)) == 36
        @test (2x^2*y)([3,2]) == 36
        @test (2x^2)(3) == 18
    end
end
