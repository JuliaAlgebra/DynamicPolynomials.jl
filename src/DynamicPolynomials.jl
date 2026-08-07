module DynamicPolynomials

using Reexport
@reexport using MultivariatePolynomials
import MultivariatePolynomials as MP

import MutableArithmetics as MA
import StarAlgebras as SA

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
<<<<<<< Updated upstream
MP.monomial_type(::Type{<:PolyType{V,M}}) where {V,M} = Monomial{V,M}
MP.monomial_type(::PolyType{V,M}) where {V,M} = Monomial{V,M}
MP.ordering(p::PolyType) = MP.ordering(MP.variable_union_type(p))
#function MP.constant_monomial(::Type{Monomial{V,M}}, vars=Variable{V,M}[]) where {V,M}
#    return Monomial{V,M}(vars, zeros(Int, length(vars)))
#end
=======

MP.constant_monomial(::Type{DPMonomial{V,M}}) where {V,M} = MP.Polynomial(
    MP.Variables{MP.Monomial}(Variable{V,M}[]),
    Int[],
)
# constant_monomial for a Variable instance
function MP.constant_monomial(v::Variable{V,M}) where {V,M}
    return MP.Polynomial(
        MP.Variables{MP.Monomial}([v]),
        [0],
    )
end
function MP.constant_monomial(::Type{Variable{V,M}}) where {V,M}
    return MP.Polynomial(
        MP.Variables{MP.Monomial}(Variable{V,M}[]),
        Int[],
    )
end
MP.monomial_type(::Type{<:DPMonomial{V,M}}) where {V,M} = DPMonomial{V,M}
MP.monomial_type(::DPMonomial{V,M}) where {V,M} = DPMonomial{V,M}
MP.monomial_type(::Type{<:Variable{V,M}}) where {V,M} = DPMonomial{V,M}
MP.monomial_type(::Variable{V,M}) where {V,M} = DPMonomial{V,M}
# MP.ordering for Variable is in var.jl

# Compute fully-parametrized SA.Term{T,A,I} type for DP variables/monomials
# We manually construct the full basis type to avoid going through
# MA.promote_operation(variables, ...) which would recurse.
function _dp_full_basis_type(::Type{Variable{V,M}}) where {V,M}
    Vars = Vector{Variable{V,M}}
    E = Vector{Int}
    P = MP.Polynomial{MP.Monomial,Vars,E}
    O = M  # monomial ordering
    return SA.MappedBasis{
        P,E,
        MP.ExponentsIterator{O,Nothing,E},
        MP.Variables{MP.Monomial,Vars},
        typeof(MP.exponents),
    }
end

function _dp_term_type(::Type{Variable{V,M}}, ::Type{T}) where {V,M,T}
    BT = _dp_full_basis_type(Variable{V,M})
    A = MA.promote_operation(MP.algebra, BT)
    I = Vector{Int}
    return SA.Term{T,A,I}
end

>>>>>>> Stashed changes
function MP.term_type(
    ::Union{TermPoly{V,M,T},Type{<:TermPoly{V,M,T}}},
) where {V,M,T}
    return _Term{V,M,T}
end
function MP.term_type(
    ::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{T},
) where {V,M,T}
<<<<<<< Updated upstream
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
=======
    return _dp_term_type(Variable{V,M}, T)
end

# term_type for Polynomial{Monomial,...} basis element (same as Variable)
function MP.term_type(::Type{<:DPMonomial{V,M}}, ::Type{T}) where {V,M,T}
    return _dp_term_type(Variable{V,M}, T)
end

# 1-arg term_type: default coefficient type Int
MP.term_type(::Type{<:DPMonomial{V,M}}) where {V,M} = _dp_term_type(Variable{V,M}, Int)
MP.term_type(::Type{<:Variable{V,M}}) where {V,M} = _dp_term_type(Variable{V,M}, Int)

MP.variables(p::AbstractArray{<:Variable}) = mergevars(MP.variables.(p))[1]
function MP.nvariables(p::Union{Variable,AbstractArray{<:Variable}})
>>>>>>> Stashed changes
    return length(MP.variables(p))
end
function MP.similar_variable(
    P::Union{PolyType{V,M},Type{<:PolyType{V,M}}},
    ::Type{Val{S}},
) where {V,M,S}
    return MP.similar_variable(P, S)
end
function MP.similar_variable(p::PolyType{V,M}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, isreal(p) ? REAL : COMPLEX)
end
function MP.similar_variable(::Type{<:PolyType{V,M}}, s::Symbol) where {V,M}
    return Variable(string(s), V, M, REAL) # we cannot infer this from the type
end

include("promote.jl")

include("operators.jl")
include("comp.jl")

include("anti_diff.jl")
include("diff.jl")
include("subs.jl")

include("div.jl")

end # module
