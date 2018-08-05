using DynamicPolynomials
using MultivariatePolynomials
using Compat
using Compat.Test
using Compat.LinearAlgebra

include("mono.jl")
include("poly.jl")
include("comp.jl")

module newmodule
    using Compat
    using Compat.Test
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
