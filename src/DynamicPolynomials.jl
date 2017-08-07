__precompile__()

module DynamicPolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one, *, +, -

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
MP.constantmonomial(p::PolyType) = Monomial(_vars(p), zeros(Int, nvars(p)))
MP.termtype(::Type{<:DMonomialLike{C}}) where {C} = Term{C, Int}
MP.termtype(::Type{<:TermPoly{C, T}}) where {C, T} = Term{C, T}
MP.termtype(::Type{<:PolyType{C}}, ::Type{T}) where {C, T} = Term{C, T}
MP.polynomial(p::PolyType) = Polynomial(p)
MP.polynomial(p::PolyType{C}, ::Type{T}) where {C, T} = Polynomial{C, T}(p)
MP.polynomialtype(::Type{<:DMonomialLike{C}}) where {C} = Polynomial{C, Int}
MP.polynomialtype(::Type{Term{C, T}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Type{T}, ::Type{<:DMonomialLike{C}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Type{<:PolyType{C}}, ::Type{T}) where {C, T} = Polynomial{C, T}
MP.vars(p::Union{PolyType, MonomialVector}) = _vars(p) # tuple(_vars(p))
MP.nvars(p::Union{PolyType, MonomialVector}) = length(_vars(p))
include("promote.jl")

include("operators.jl")
include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")

include("show.jl")

end # module
