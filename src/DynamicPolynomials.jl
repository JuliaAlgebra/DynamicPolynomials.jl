__precompile__()

module DynamicPolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one, *, +, -

using MultivariatePolynomials
const MP = MultivariatePolynomials

#const PolyType{C} = Union{DMonomialLike{C}, RationalPoly{C}}
#iscomm{C}(::PolyType{C}) = C
#zero{C}(p::PolyType{C}) = zero(typeof(p))
#one{C}(p::PolyType{C}) = one(typeof(p))

include("var.jl")
include("mono.jl")
const DMonomialLike{C} = Union{Monomial{C}, PolyVar{C}}
include("term.jl")
include("monovec.jl")
include("poly.jl")
const TermPoly{C, T} = Union{Term{C, T}, Polynomial{C, T}}
const PolyType{C} = Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}
MP.termtype{C}(::Type{<:DMonomialLike{C}}) = Term{C, Int}
MP.termtype{TT<:Term}(::Type{TT}) = TT
MP.termtype{C, T}(::Type{Polynomial{C, T}}) = Term{C, T}
MP.term(α, p::PolyType) = α * Monomial(_vars(p), zeros(Int, nvars(p)))
MP.polynomial(p::PolyType) = Polynomial(p)
MP.polynomial{C, T}(p::PolyType{C}, ::Type{T}) = Polynomial{C, T}(p)
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
