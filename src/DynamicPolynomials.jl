module DynamicPolynomials

import Future # For `copy!`

using Reexport
@reexport using MultivariatePolynomials
import MultivariatePolynomials as MP

import MutableArithmetics as MA

include("var.jl")
#const CommutativeVariable{O,M} = Variable{Commutative{O},M}
#const NonCommutativeVariable{O,M} = Variable{NonCommutative{O},M}
include("mono.jl")
const DMonomialLike{V,M} = Union{Monomial{V,M},Variable{V,M}}
MA.mutability(::Type{<:Monomial{<:Commutative}}) = MA.IsMutable()
MA.mutability(::Type{<:Monomial{<:NonCommutative}}) = MA.IsNotMutable()
const _Term{V,M,T} = MP.Term{T,Monomial{V,M}}
function __add_variables!(t::_Term, allvars, map)
    return __add_variables!(MP.monomial(t), allvars, map)
end
include("monomial_vector.jl")
include("poly.jl")
MA.mutability(::Type{<:Polynomial}) = MA.IsMutable()
const TermPoly{V,M,T} = Union{_Term{V,M,T},Polynomial{V,M,T}}
const PolyType{V,M} =
    Union{Polynomial{V,M},_Term{V,M},Monomial{V,M},Variable{V,M}}
function MP.variable_union_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
) where {V,M}
    return Variable{V,M}
end
MP.constant_monomial(::Type{<:PolyType{V,M}}) where {V,M} = Monomial{V,M}()
function MP.constant_monomial(p::PolyType)
    return Monomial(copy(MP.variables(p)), zeros(Int, nvariables(p)))
end
MP.monomial_type(::Type{<:PolyType{V,M}}) where {V,M} = Monomial{V,M}
MP.monomial_type(::PolyType{V,M}) where {V,M} = Monomial{V,M}
#function MP.constant_monomial(::Type{Monomial{V,M}}, vars=Variable{V,M}[]) where {V,M}
#    return Monomial{V,M}(vars, zeros(Int, length(vars)))
#end
function MP.term_type(
    ::Union{TermPoly{V,M,T},Type{<:TermPoly{V,M,T}}},
) where {V,M,T}
    return _Term{V,M,T}
end
function MP.term_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{T},
) where {V,M,T}
    return _Term{V,M,T}
end
MP.term_type(::Type{Polynomial{V,M}}) where {V,M} = _Term{V,M}
MP.polynomial_type(::Type{_Term{V,M}}) where {V,M} = Polynomial{V,M}
MP.polynomial_type(::Type{_Term{V,M,T}}) where {T,V,M} = Polynomial{V,M,T}
function MP.polynomial_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{T},
) where {V,M,T}
    return Polynomial{V,M,T}
end
MP.variables(p::AbstractArray{<:PolyType}) = mergevars(MP.variables.(p))[1]
function MP.nvariables(
    p::Union{PolyType,MonomialVector,AbstractArray{<:PolyType}},
)
    return length(MP.variables(p))
end
function MP.similar_variable(
    P::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{Val{S}},
) where {V,M,S}
    return MP.similar_variable(P, S)
end
function MP.similar_variable(p::PolyType{V,M}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, iscomplex(p) ? cpFull : cpNone)
end
function MP.similar_variable(::Type{<:PolyType{V,M}}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, cpNone) # we cannot infer this from the type,
end

include("promote.jl")

include("operators.jl")
include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")

end # module
