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

@testset "Issue #79 and Issue #80" begin
    @polyvar x[1:2]
    p1 = x[1] * 0.0 + x[2] * 0
    p2 = ( x[1] + x[2] ) * 0.0
    @test variables(p1) == x 
    @test variables(p1) == variables(p2)

    @polyvar χ[1:4]
    poly_array_int = [1 0 1 0; 0 1 0 1]*χ;
    poly_array_float = Float64.([1 0 1 0; 0 1 0 1])*χ;

    vars_array_int = variables.(poly_array_int);
    vars_array_float = variables.(poly_array_float)
    @test isa(vars_array_int, Vector{<:Vector{<:PolyVar}})
    @test isa(vars_array_float, Vector{<:Vector{<:PolyVar}})
    @test vars_array_int == vars_array_float
    @test vars_array_int[1] == vars_array_float[2] == χ
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
