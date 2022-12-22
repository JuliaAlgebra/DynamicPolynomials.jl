module DynamicPolynomials

import Future # For `copy!`

using Reexport
@reexport using MultivariatePolynomials
const MP = MultivariatePolynomials

import MutableArithmetics
const MA = MutableArithmetics

using DataStructures

include("var.jl")
#const CommutativeVariable{O,M} = Variable{Commutative{O},M}
#const NonCommutativeVariable{O,M} = Variable{NonCommutative{O},M}
include("mono.jl")
const DMonomialLike{V,M} = Union{Monomial{V,M},Variable{V,M}}
MA.mutability(::Type{<:Monomial}) = MA.IsMutable()
const _Term{V,M,T} = MP.Term{T,Monomial{V,M}}
include("monomial_vector.jl")
include("poly.jl")
MA.mutability(::Type{<:Polynomial}) = MA.IsMutable()
const TermPoly{V,M,T} = Union{_Term{V,M,T},Polynomial{V,M,T}}
const PolyType{V,M} = Union{Polynomial{V,M},_Term{V,M},Monomial{V,M},Variable{V,M}}
function MP.variable_union_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
) where {V,M}
    return Variable{V,M}
end
MP.constant_monomial(::Type{<:PolyType{V,M}}) where {V,M} = Monomial{V,M}()
function MP.constant_monomial(p::PolyType)
    return Monomial(_vars(p), zeros(Int, nvariables(p)))
end
MP.monomial_type(::Type{<:PolyType{V,M}}) where {V,M} = Monomial{V,M}
MP.monomial_type(::PolyType{V,M}) where {V,M} = Monomial{V,M}
#function MP.constant_monomial(::Type{Monomial{V,M}}, vars=Variable{V,M}[]) where {V,M}
#    return Monomial{V,M}(vars, zeros(Int, length(vars)))
#end
function MP.term_type(::Union{TermPoly{V,M,T},Type{<:TermPoly{V,M,T}}}) where {V,M,T}
    return Term{V,M,T}
end
function MP.term_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{T},
) where {V,M,T}
    return Term{V,M,T}
end
MP.term_type(::Type{Polynomial{V,M}}) where {V,M} = Term{V,M}
MP.polynomial_type(::Type{_Term{V,M}}) where {V,M} = Polynomial{V,M}
MP.polynomial_type(::Type{_Term{V,M,T}}) where {T,V,M} = Polynomial{V,M,T}
function MP.polynomial_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{T},
) where {V,M,T}
    return Polynomial{V,M,T}
end
_vars(p::AbstractArray{<:PolyType}) = mergevars(_vars.(p))[1]
function MP.variables(
    p::Union{PolyType,MonomialVector,AbstractArray{<:PolyType}},
)
    return _vars(p)
end # tuple(_vars(p))
function MP.nvariables(
    p::Union{PolyType,MonomialVector,AbstractArray{<:PolyType}},
)
    return length(_vars(p))
end
function MP.similar_variable(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{Val{S}},
) where {V,M,S}
    return Variable{V,M}(string(S))
end
function MP.similar_variable(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    s::Symbol,
) where {V,M}
    return Variable{V,M}(string(s))
end

include("promote.jl")

include("operators.jl")
include("comp.jl")
#
#include("diff.jl")
#include("subs.jl")
#
#include("div.jl")

end # module
