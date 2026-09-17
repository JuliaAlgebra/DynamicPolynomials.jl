module DynamicPolynomials

using Reexport
@reexport using MultivariatePolynomials
import MultivariatePolynomials as MP

import MutableArithmetics as MA
import StarAlgebras as SA

include("var.jl")

# Monomials are now MP.Polynomial{MP.Monomial, V, E}
# No separate Monomial struct — variables and exponents live in the basis element.

# Convenience alias for the specific monomial type with DP variables
const DPMonomial{V,M} = MP.Polynomial{MP.Monomial,Vector{Variable{V,M}},Vector{Int}}

function MP.variable_union_type(
    ::Union{Variable{V,M},Type{<:Variable{V,M}}},
) where {V,M}
    return Variable{V,M}
end
function MP.variable_union_type(
    ::Union{DPMonomial{V,M},Type{<:DPMonomial{V,M}}},
) where {V,M}
    return Variable{V,M}
end

MP.constant_monomial(::Type{DPMonomial{V,M}}) where {V,M} = MP.Polynomial(
    MP.Variables{MP.Monomial}(Variable{V,M}[]),
    Int[],
)
MP.monomial_type(::Type{<:Variable{V,M}}) where {V,M} = DPMonomial{V,M}
MP.monomial_type(::Variable{V,M}) where {V,M} = DPMonomial{V,M}
# MP.ordering for Variable is in var.jl

MP.variables(p::AbstractArray{<:Variable}) = mergevars(MP.variables.(p))[1]
function MP.nvariables(p::Union{Variable,AbstractArray{<:Variable}})
    return length(MP.variables(p))
end
function MP.similar_variable(
    P::Union{Variable{V,M},Type{<:Variable{V,M}}},
    ::Type{Val{S}},
) where {V,M,S}
    return MP.similar_variable(P, S)
end
function MP.similar_variable(p::Variable{V,M}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, isreal(p) ? REAL : COMPLEX)
end
function MP.similar_variable(::Type{<:Variable{V,M}}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, REAL)
end

# Create monomial from variable: Variable → Polynomial{Monomial,...}
function Base.convert(::Type{DPMonomial{V,M}}, x::Variable{V,M}) where {V,M}
    return MP.Polynomial{MP.Monomial}(x)
end

# monomial(vars, exps) constructs a Polynomial{Monomial,...}
function MP.monomial(vars::Vector{Variable{V,M}}, z::Vector{Int}) where {V,M}
    @assert !iscomm(V) || issorted(vars, rev = true)
    return MP.Polynomial(MP.Variables{MP.Monomial}(vars), z)
end

# exponents for variables is in var.jl

include("comp.jl")
include("promote.jl")

# Variable power → monomial
Base.:(^)(x::Variable{V,M}, i::Int) where {V,M} = MP.Polynomial(
    MP.Variables{MP.Monomial}([x]),
    [i],
)

# Variable + Variable → uses term + term → AlgebraElement
Base.:(+)(x::Variable, y::Variable) = MP.term(x) + MP.term(y)
Base.:(-)(x::Variable, y::Variable) = MP.term(x) - MP.term(y)

# Short names for the default commutative algebra and coefficient storage.
const _DefaultVariable =
    Variable{Commutative{CreationOrder},MP.Graded{MP.LexOrder}}
const _DefaultAlgebra =
    typeof(MP.algebra(MP.FullBasis{MP.Monomial}(_DefaultVariable[])))
const Polynomial{T} = SA.AlgebraElement{
    T,
    _DefaultAlgebra,
    SA.SparseCoefficients{
        Vector{Int},
        T,
        Vector{Vector{Int}},
        Vector{T},
        MP.Graded{MP.LexOrder},
    },
}
const Term{T} = SA.Term{T,_DefaultAlgebra,Vector{Int}}

# Julia searches for type aliases only in the outer type's defining module.
function Base.show(io::IO, ::Type{Polynomial{T}}) where {T}
    print(io, "DynamicPolynomials.Polynomial{")
    show(io, T)
    return print(io, "}")
end
function Base.show(io::IO, ::Type{Term{T}}) where {T}
    print(io, "DynamicPolynomials.Term{")
    show(io, T)
    return print(io, "}")
end

end # module
