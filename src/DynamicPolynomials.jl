module DynamicPolynomials

import Future # For `copy!`

using Reexport
@reexport using MultivariatePolynomials
const MP = MultivariatePolynomials

import MutableArithmetics
const MA = MutableArithmetics

using DataStructures

include("var.jl")
include("mono.jl")
const DMonomialLike{C} = Union{Monomial{C}, PolyVar{C}}
MA.mutability(::Type{<:Monomial}) = MA.IsMutable()
include("term.jl")
MA.mutability(::Type{Term{C, T}}) where {C, T} = MA.mutability(T)
include("monovec.jl")
include("poly.jl")
MA.mutability(::Type{<:Polynomial}) = MA.IsMutable()
const TermPoly{C, T} = Union{Term{C, T}, Polynomial{C, T}}
const PolyType{C} = Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}
MP.variable_union_type(::Union{PolyType{C}, Type{<:PolyType{C}}}) where {C} = PolyVar{C}
MP.constantmonomial(::Type{<:PolyType{C}}) where {C} = Monomial{C}()
MP.constantmonomial(p::PolyType) = Monomial(_vars(p), zeros(Int, nvariables(p)))
MP.monomialtype(::Type{<:PolyType{C}}) where C = Monomial{C}
MP.monomialtype(::PolyType{C}) where C = Monomial{C}
#function MP.constantmonomial(::Type{Monomial{C}}, vars=PolyVar{C}[]) where {C}
#    return Monomial{C}(vars, zeros(Int, length(vars)))
#end
MP.termtype(::Union{TermPoly{C, T}, Type{<:TermPoly{C, T}}}) where {C, T} = Term{C, T}
MP.termtype(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Term{C, T}
MP.termtype(::Type{Polynomial{C}}) where {C} = Term{C}
MP.polynomialtype(::Type{Term{C}}) where {C} = Polynomial{C}
MP.polynomialtype(::Type{Term{C, T}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Type{T}, ::Type{<:DMonomialLike{C}}) where {T, C} = Polynomial{C, T}
MP.polynomialtype(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Polynomial{C, T}
_vars(p::AbstractArray{<:PolyType}) = mergevars(_vars.(p))[1]
MP.variables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = _vars(p) # tuple(_vars(p))
MP.nvariables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = length(_vars(p))
MP.similarvariable(p::Type{<:PolyType{C}}, ::Type{Val{V}}) where {C, V} = PolyVar{C}(string(V))
MP.similarvariable(p::Type{<:PolyType{C}}, V::Symbol) where {C} = PolyVar{C}(string(V))
MP.similarvariable(p::PolyVar{C}, ::Type{Val{V}}) where {C,V} = PolyVar{C}(string(V), p.kind)
MP.similarvariable(p::PolyVar{C}, V::Symbol) where {C} = similarvariable(p, Val{V})

include("promote.jl")

include("operators.jl")
include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")

end # module
