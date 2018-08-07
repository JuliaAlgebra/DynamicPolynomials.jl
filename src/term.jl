export Term

struct Term{C, T} <: AbstractTerm{T}
    α::T
    x::Monomial{C}
end

iscomm(::Type{Term{C, T}}) where {C, T} = C

Term{C, T}(t::Term{C, T}) where {C, T} = t
Term{C, T}(t::Term{C}) where {C, T} = Term{C, T}(T(t.α), t.x)
Term(t::Term) = t

Term{C, T}(x::Monomial{C}) where {C, T} = Term{C, T}(one(T), x)
Term{C}(x::Monomial{C}) where C = Term{C, Int}(x)
Term(x::Monomial{C}) where C = Term{C}(x)

Term{C, T}(x::PolyVar{C}) where {C, T} = Term{C, T}(Monomial{C}(x))
Term{C}(x::PolyVar{C}) where C = Term{C}(Monomial{C}(x))
Term(x::PolyVar{C}) where C = Term{C}(x)

Term{C, T}(α) where {C, T} = Term{C}(T(α))
Term{C}(α::T) where {C, T} = Term{C, T}(α, Monomial{C}())

Base.broadcastable(t::Term) = Ref(t)
#(::Type{TermContainer{C}}){C}(x::PolyVar{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(x::Monomial{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(t::TermContainer{C}) = t
#TermContainer(x::PolyVar) = Term(x)
#TermContainer(x::Monomial) = Term(x)
#TermContainer(t::TermContainer) = t

#Base.convert{C, T}(::Type{TermContainer{C, T}}, x::Union{Monomial{C},PolyVar{C}}) = Term{C, T}(x)
#Base.convert{C, T}(::Type{TermContainer{C, T}}, α::T) = Term{C, T}(α, Monomial{C}())
#Base.convert{C, S, T}(::Type{TermContainer{C, T}}, α::S) = TermContainer{C, T}(T(α))
#(::Type{TermContainer{C}}){C, T}(α::T) = TermContainer{C, T}(α)

#Base.convert{C, T}(::Type{TermContainer{C, T}}, t::Term{C}) = Term{C, T}(t)

Base.copy(t::T) where {T<:Term} = T(copy(t.α), copy(t.x))

MP.coefficient(t::Term) = t.α
MP.monomial(t::Term) = t.x
_vars(t) = _vars(t.x)
