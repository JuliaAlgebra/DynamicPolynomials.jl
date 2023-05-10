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
include("monomial_vector.jl")
include("poly.jl")
MA.mutability(::Type{<:Polynomial}) = MA.IsMutable()
const TermPoly{C, T} = Union{Term{C, T}, Polynomial{C, T}}
const PolyType{C} = Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}
MP.variable_union_type(::Union{PolyType{C}, Type{<:PolyType{C}}}) where {C} = PolyVar{C}
MP.constant_monomial(::Type{<:PolyType{C}}) where {C} = Monomial{C}()
MP.constant_monomial(p::PolyType) = Monomial(_vars(p), zeros(Int, nvariables(p)))
MP.monomial_type(::Type{<:PolyType{C}}) where C = Monomial{C}
MP.monomial_type(::PolyType{C}) where C = Monomial{C}
#function MP.constant_monomial(::Type{Monomial{C}}, vars=PolyVar{C}[]) where {C}
#    return Monomial{C}(vars, zeros(Int, length(vars)))
#end
MP.term_type(::Union{TermPoly{C, T}, Type{<:TermPoly{C, T}}}) where {C, T} = Term{C, T}
MP.term_type(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Term{C, T}
MP.term_type(::Type{Polynomial{C}}) where {C} = Term{C}
MP.polynomial_type(::Type{Term{C}}) where {C} = Polynomial{C}
MP.polynomial_type(::Type{Term{C, T}}) where {T, C} = Polynomial{C, T}
MP.polynomial_type(::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{T}) where {C, T} = Polynomial{C, T}
_vars(p::AbstractArray{<:PolyType}) = mergevars(_vars.(p))[1]
MP.variables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = _vars(p) # tuple(_vars(p))
MP.nvariables(p::Union{PolyType, MonomialVector, AbstractArray{<:PolyType}}) = length(_vars(p))
MP.similar_variable(p::Union{PolyType{C}, Type{<:PolyType{C}}}, ::Type{Val{V}}) where {C, V} = PolyVar{C}(string(V))
MP.similar_variable(p::Union{PolyType{C}, Type{<:PolyType{C}}}, V::Symbol) where {C} = PolyVar{C}(string(V))

include("promote.jl")

include("operators.jl")
include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")

end # module
