using DynamicPolynomials
using MultivariatePolynomials
using Test
using LinearAlgebra

function alloc_test_lt(f, n)
    f() # compile
    @test n >= @allocated f()
end

# TODO move to MP
@testset "See https://github.com/jump-dev/SumOfSquares.jl/issues/388" begin
    @polyvar x[1:3]
    p = sum(x)
    v = map(_ -> 1, x)
    # I get 208 but let's give some margin
    alloc_test_lt(() -> substitute(Eval(), p, x => v), 300)
    alloc_test_lt(() -> p(x => v), 300)
    alloc_test_lt(() -> substitute(Eval(), p, v), 0)
    alloc_test_lt(() -> p(v), 0)
    err = ErrorException("Cannot evaluate a polynomial of `3` variables with only `2` values.")
    @test_throws err p([1, 2])
end

@testset "Issue #70" begin
    @ncpolyvar y0 y1 x0 x1
    p = x1 * x0 * x1
    @test subs(p, x0 => y0, x1 => y1) == y1 * y0 * y1
    @test subs(p, x0 => 1) == x1^2
    @test p(x0 => y0, x1 => y1) == y1 * y0 * y1
end

# https://github.com/JuliaAlgebra/DynamicPolynomials.jl/issues/141
@testset "Issue #141" begin
    @polyvar x
    m = x^2
    q = m^2
    @test variables(q) !== variables(m)
end

@testset "Issue #79, Issue #80 and Issue #92" begin
    @polyvar x[1:2]
    p1 = x[1] * 0.0 + x[2] * 0
    p2 = (x[1] + x[2]) * 0.0
    @test variables(p1) == x
    @test variables(p1) == variables(p2)

    @polyvar χ[1:4]
    poly_array_int = [1 0 1 0; 0 1 0 1] * χ
    poly_array_float = Float64.([1 0 1 0; 0 1 0 1]) * χ

    vars_array_int = variables.(poly_array_int)
    vars_array_float = variables.(poly_array_float)
    @test isa(vars_array_int, Vector{<:Vector{<:Variable}})
    @test isa(vars_array_float, Vector{<:Vector{<:Variable}})
    @test vars_array_int == vars_array_float
    @test vars_array_int[1] == vars_array_float[2] == χ

    # This is for Issue #92
    p3 = subs(p1, x => x)
    @test variables(p3) == x

    p4 = x[1] + x[2]
    @test variables(subs(p4, x[2] => 1)) == [x[1]]
    @test variables(subs(p4, x[1] => 1)) == [x[2]]
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
    @test isa(x, DynamicPolynomials.Variable)
    @test isa(y, DynamicPolynomials.Variable)
end
end

@testset "Issue #166: promote_operation with Any" begin
    DynamicPolynomials.@polyvar x
    V = typeof(x)
    @test promote_type(V, Any) == Any
end

include("mvp.jl")
