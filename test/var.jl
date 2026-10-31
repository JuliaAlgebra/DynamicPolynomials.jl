module TestVariableMacros

using Test
import DynamicPolynomials as DP

const real_vars = DP.@polyvar x xs[1:10] constant = true
const nc_vars = DP.@ncpolyvar y ys[1:2] constant = true
const complex_vars = DP.@complex_polyvar z zs[1:2, 1:3] constant = true
DP.@polyvar ordinary constant = false

@testset "Constant variable bindings" begin
    for name in (:x, :xs, :y, :ys, :z, :zs)
        @test isconst(@__MODULE__, name)
    end
    @test !isconst(@__MODULE__, :ordinary)
    @test real_vars === (x, xs)
    @test nc_vars === (y, ys)
    @test complex_vars === (z, zs)
    @test DP.name.(xs) == ["xs[$i]" for i in 1:10]
    @test size(zs) == (2, 3)
    @test DP.is_commutative(typeof(x))
    @test !DP.is_commutative(typeof(y))
    @test isreal(x)
    @test !isreal(z)
    @test all(!isreal, zs)
    @test_throws ErrorException macroexpand(
        @__MODULE__,
        :(DP.@polyvar invalid constant = 1),
    )
end

end
