using DynamicPolynomials
using MultivariatePolynomials
using Test
using LinearAlgebra

include("mono.jl")
include("poly.jl")
include("comp.jl")

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
