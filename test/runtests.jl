using DynamicPolynomials
using MultivariatePolynomials
using Test
using LinearAlgebra

# TODO move to MP
@testset "Noncommutative quadratic" begin
    @ncpolyvar x[1:2]
    Q = Hermitian([1 2 + 3im; 2 - 3im 4])
    p = 1x[1]^2 + (2 + 3im) * x[1] * x[2] + (2 - 3im) * x[2] * x[1] + 4x[2]^2
    @test polynomial(Q, x) == p
    @test polynomial(Q, monovec(x)) == p
end

include("mono.jl")
include("poly.jl")
include("comp.jl")
include("mutable_arithmetics.jl")

# TODO move to MultivariatePolynomials.jl
@testset "Subs with no variables" begin
    @polyvar x
    t = convert(termtype(x, Int), 3)
    @test t == @inferred subs(t, x => x + 1)
    @test t == @inferred subs(t, x => x + 1.0)
    @test t == @inferred subs(t, x => 1x)
    @test t == @inferred subs(t, x => 1.0x)
    @test t == @inferred subs(t, x => 1.0)
end

module newmodule
    using Test
    import DynamicPolynomials
    @testset "Polyvar macro hygiene" begin
        # Verify that the @polyvar macro works when the package has been activated
        # with `import` instead of `using`.
        DynamicPolynomials.@polyvar x y
        @test isa(x, DynamicPolynomials.PolyVar)
        @test isa(y, DynamicPolynomials.PolyVar)
    end
end

include("mvp.jl")
