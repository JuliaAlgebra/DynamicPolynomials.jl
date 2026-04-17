module DynamicPolynomials

using Reexport
@reexport using MultivariatePolynomials
import MultivariatePolynomials as MP

import MutableArithmetics as MA
import StarAlgebras as SA

include("var.jl")
include("mono.jl")
const DMonomialLike{V,M} = Union{Monomial{V,M},Variable{V,M}}
MA.mutability(::Type{<:Monomial{<:Commutative}}) = MA.IsMutable()
MA.mutability(::Type{<:Monomial{<:NonCommutative}}) = MA.IsNotMutable()

# Term type comes from MP (SA.Term)
const _Term{V,M,T} = MP.Term{T,Monomial{V,M}}

include("monomial_vector.jl")

function MP.variable_union_type(
    ::Union{DMonomialLike{V,M},Type{<:DMonomialLike{V,M}}},
) where {V,M}
    return Variable{V,M}
end
function MP.variable_union_type(
    ::Union{_Term{V,M},Type{<:_Term{V,M}}},
) where {V,M}
    return Variable{V,M}
end
MP.constant_monomial(::Type{<:DMonomialLike{V,M}}) where {V,M} = Monomial{V,M}()
function MP.constant_monomial(p::DMonomialLike)
    return Monomial(copy(MP.variables(p)), zeros(Int, nvariables(p)))
end
MP.monomial_type(::Type{<:DMonomialLike{V,M}}) where {V,M} = Monomial{V,M}
MP.monomial_type(::DMonomialLike{V,M}) where {V,M} = Monomial{V,M}
MP.monomial_type(::Type{<:_Term{V,M}}) where {V,M} = Monomial{V,M}
MP.monomial_type(::_Term{V,M}) where {V,M} = Monomial{V,M}
MP.ordering(p::DMonomialLike) = MP.ordering(MP.variable_union_type(p))

function MP.term_type(
    ::Union{DMonomialLike{V,M},Type{<:DMonomialLike{V,M}}},
    ::Type{T},
) where {V,M,T}
    return _Term{V,M,T}
end

# polynomial_type returns AlgebraElement type
# We compute the FullBasis type directly to avoid recursion through _polynomial_type
function _dp_full_basis_type(::Type{Monomial{V,M}}) where {V,M}
    Vars = Vector{Variable{V,M}}
    E = Vector{Int}
    P = MP.Polynomial{MP.Monomial,Vars,E}
    O = MP.ordering(Variable{V,M})
    return SA.MappedBasis{
        P, E,
        MP.ExponentsIterator{O,Nothing,E},
        MP.Variables{MP.Monomial,Vars},
        typeof(MP.exponents),
    }
end

function MP.polynomial_type(
    ::Union{DMonomialLike{V,M},Type{<:DMonomialLike{V,M}}},
    ::Type{T},
) where {V,M,T}
    return MP.algebra_element_type(Vector{T}, _dp_full_basis_type(Monomial{V,M}))
end
function MP.polynomial_type(::Type{_Term{V,M,T}}) where {V,M,T}
    return MP.algebra_element_type(Vector{T}, _dp_full_basis_type(Monomial{V,M}))
end
function MP.polynomial_type(::Type{_Term{V,M}}) where {V,M}
    return MP.algebra_element_type(Vector{Any}, _dp_full_basis_type(Monomial{V,M}))
end

MP.variables(p::AbstractArray{<:DMonomialLike}) = mergevars(MP.variables.(p))[1]
function MP.nvariables(
    p::Union{DMonomialLike,MonomialVector,AbstractArray{<:DMonomialLike}},
)
    return length(MP.variables(p))
end
function MP.similar_variable(
    P::Union{DMonomialLike{V,M},Type{<:DMonomialLike{V,M}}},
    ::Type{Val{S}},
) where {V,M,S}
    return MP.similar_variable(P, S)
end
function MP.similar_variable(p::DMonomialLike{V,M}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, isreal(p) ? REAL : COMPLEX)
end
function MP.similar_variable(::Type{<:DMonomialLike{V,M}}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, REAL)
end

include("promote.jl")

# Monomial multiplication (cmult.jl and ncmult.jl)
include("mult.jl")

include("comp.jl")

include("diff.jl")
include("subs.jl")

include("div.jl")

end # module
