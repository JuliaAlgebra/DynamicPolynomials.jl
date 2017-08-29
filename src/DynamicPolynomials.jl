__precompile__()

module DynamicPolynomials

import Base: length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one, *, +, -

using MultivariatePolynomials
const MP = MultivariatePolynomials

#const PolyType{C} = Union{DMonomialLike{C}, RationalPoly{C}}
#iscomm(::PolyType{C}) where {C} = C
#zero(p::PolyType{C}) where {C} = zero(typeof(p))
#one(p::PolyType{C}) where {C} = one(typeof(p))

include("var.jl")
include("mono.jl")
const DMonomialLike{C} = Union{Monomial{C}, PolyVar{C}}
include("term.jl")
include("monovec.jl")
include("poly.jl")
const TermPoly{C, T} = Union{Term{C, T}, Polynomial{C, T}}
const PolyType{C} = Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}
MP.constantmonomial(::Type{<:PolyType{C}}) where {C} = Monomial{C}()
MP.constantmonomial(p::PolyType) = Monomial(_vars(p), zeros(Int, nvariables(p)))
MP.monomialtype(::Type{<:PolyType{C}}) where C = Monomial{C}
MP.monomialtype(::PolyType{C}) where C = Monomial{C}
MP.termtype(::Union{TermPoly{C, T}, Type{<:TermPoly{C, T}}}) where {C, T} = Term{C, T}
MP.termtype(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Term{C, T}
MP.polynomial(p::PolyType) = Polynomial(p)
MP.polynomial(p::PolyType{C}, ::Type{T}) where {C, T} = Polynomial{C, T}(p)
MP.polynomialtype(::Type{Term{C, T}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Type{T}, ::Type{<:DMonomialLike{C}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Polynomial{C, T}
_vars(p::AbstractArray{<:PolyType}) = mergevars(_vars.(p))[1]
MP.variables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = _vars(p) # tuple(_vars(p))
MP.nvariables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = length(_vars(p))
MP.similarvariable(p::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{Val{V}}) where {C, V} = PolyVar{C}(string(V))
include("promote.jl")

include("operators.jl")
include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")


end # module
