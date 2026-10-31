using Test
import DynamicPolynomials as DP
import MultivariatePolynomials as MP
import StarAlgebras as SA

@testset "Compact polynomial type display" begin
    DP.@polyvar x y
    for T in (Int, Float64, ComplexF64)
        term = MP.term(one(T), x^2)
        polynomial = one(T) * x + one(T) * y
        coefficients = SA.SparseCoefficients(
            (term.index,),
            (MP.coefficient(term),),
            MP.ordering(term),
        )
        single = SA.AlgebraElement(coefficients, parent(term))
        @test typeof(polynomial) === DP.Polynomial{T}
        @test typeof(term) === DP.Term{T}
        @test typeof(single) === DP.SingleTermPolynomial{T}
        for (value, name) in (
            (polynomial, "Polynomial"),
            (term, "Term"),
            (single, "SingleTermPolynomial"),
        )
            expected = "DynamicPolynomials.$name{$T}"
            @test sprint(show, typeof(value)) == expected
            @test sprint(show, Tuple{typeof(value)}) == "Tuple{$expected}"
            @test sprint(show, Vector{typeof(value)}) == "Vector{$expected}"
        end
    end
    for (family, name) in (
        (DP.Polynomial, "Polynomial"),
        (DP.Term, "Term"),
        (DP.SingleTermPolynomial, "SingleTermPolynomial"),
    )
        @test sprint(show, family) == "DynamicPolynomials.$name"
        @test sprint(show, family{Float64}) ==
              "DynamicPolynomials.$name{Float64}"
        bounded = family{T} where {T<:Real}
        @test sprint(show, bounded) == "DynamicPolynomials.$name{<:Real}"
    end
    first_type = DP.Term{Int}
    second_type = DP.Term{Float64}
    union_display = sprint(show, Union{first_type,second_type})
    @test startswith(union_display, "Union{")
    @test occursin(sprint(show, first_type), union_display)
    @test occursin(sprint(show, second_type), union_display)
    @test sprint(show, Union{}) == "Union{}"
    @test isempty(filter(m -> m.name == :show, Test.detect_unbound_args(DP)))
end
