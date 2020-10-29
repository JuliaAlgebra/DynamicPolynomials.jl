using DynamicPolynomials
using MultivariatePolynomials
using Test
using LinearAlgebra

# TODO move to MP
@testset "Issue #70" begin
    @ncpolyvar y0 y1 x0 x1
    p = x1 * x0 * x1
    @test subs(p, x0 => y0, x1 => y1) == y1 * y0 * y1
    @test subs(p, x0 => 1) == x1^2
    @test p(x0 => y0, x1 => y1) == y1 * y0 * y1
end

include("mono.jl")
include("poly.jl")
include("comp.jl")
include("mutable_arithmetics.jl")

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
