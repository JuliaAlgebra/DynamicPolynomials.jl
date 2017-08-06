export Term

type Term{C, T} <: AbstractTerm{T}
    α::T
    x::Monomial{C}
end

iscomm{C, T}(::Type{Term{C, T}}) = C
Term(t::Term) = t
(::Type{Term{C}}){C}(x::Monomial{C}) = Term{C, Int}(x)
(::Type{Term{C}}){C}(x::PolyVar{C}) = Term{C}(Monomial{C}(x))
Term{C}(x::Monomial{C}) = Term{C}(x)
Term{C}(x::PolyVar{C}) = Term{C}(x)
#(::Type{TermContainer{C}}){C}(x::PolyVar{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(x::Monomial{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(t::TermContainer{C}) = t
#TermContainer(x::PolyVar) = Term(x)
#TermContainer(x::Monomial) = Term(x)
#TermContainer(t::TermContainer) = t
(::Type{Term{C}}){C, T}(α::T) = Term{C, T}(α, Monomial{C}())
Base.convert{C, T}(::Type{Term{C, T}}, t::Term{C, T}) = t
Base.convert{C, T}(::Type{Term{C, T}}, t::Term{C}) = Term{C, T}(T(t.α), t.x)
Base.convert{C, T}(::Type{Term{C, T}}, x::Monomial{C}) = Term{C, T}(one(T), x)
Base.convert{C, T}(::Type{Term{C, T}}, x::PolyVar{C}) = Term{C, T}(Monomial{C}(x))
Base.convert{C, T}(::Type{Term{C, T}}, α) = Term{C}(T(α))

#Base.convert{C, T}(::Type{TermContainer{C, T}}, x::Union{Monomial{C},PolyVar{C}}) = Term{C, T}(x)
#Base.convert{C, T}(::Type{TermContainer{C, T}}, α::T) = Term{C, T}(α, Monomial{C}())
#Base.convert{C, S, T}(::Type{TermContainer{C, T}}, α::S) = TermContainer{C, T}(T(α))
#(::Type{TermContainer{C}}){C, T}(α::T) = TermContainer{C, T}(α)

#Base.convert{C, T}(::Type{TermContainer{C, T}}, t::Term{C}) = Term{C, T}(t)

Base.convert(::Type{Any}, t::Term) = t
Base.copy{T<:Term}(t::T) = T(copy(t.α), copy(t.x))

MP.coefficient(t::Term) = t.α
MP.monomial(t::Term) = t.x
_vars(t) = _vars(t.x)

#eltype(t::Term) = T
Base.length(::Term) = 1
Base.isempty(::Term) = false
Base.start(::Term) = false
Base.done(::Term, state) = state
Base.next(t::Term, state) = (t, true)
Base.getindex(t::Term, I::Int) = t
