__precompile__()

module DynamicPolynomials

import Base: show, length, getindex, vect, isless, isempty, start, done, next, convert, dot, copy, eltype, zero, one

using MultivariatePolynomials

abstract type DMonomialLike{C} end

abstract type TermType{C, T} end
abstract type TermContainer{C, T} <: TermType{C, T} end

const PolyType{C} = Union{DMonomialLike{C}, TermType{C}, RationalPoly{C}}
iscomm{C}(::PolyType{C}) = C
zero{C}(p::PolyType{C}) = zero(typeof(p))
one{C}(p::PolyType{C}) = one(typeof(p))

include("mono.jl")
include("poly.jl")
include("rational.jl")
include("measure.jl")
include("exp.jl")
include("promote.jl")

include("comp.jl")

include("alg.jl")
include("calg.jl")
include("ncalg.jl")

include("diff.jl")
include("subs.jl")
include("algebraicset.jl")
include("norm.jl")

include("div.jl")

include("show.jl")

end # module
