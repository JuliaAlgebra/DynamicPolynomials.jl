export Term

struct Term{C,T} <: AbstractTerm{T}
    α::T
    x::Monomial{C}
end
MP.term(α, mono::Monomial) = Term(α, mono)

iscomm(::Type{Term{C,T}}) where {C,T} = C

Base.convert(::Type{Term{C,T}}, t::Term{C,T}) where {C,T} = t
function Base.convert(::Type{Term{C,T}}, t::Term{C}) where {C,T}
    return Term{C,T}(T(t.α), t.x)
end
Term(t::Term) = t

function Base.convert(::Type{Term{C,T}}, x::Monomial{C}) where {C,T}
    return Term{C,T}(one(T), x)
end
Term{C}(x::Monomial{C}) where {C} = convert(Term{C,Int}, x)
Term(x::Monomial{C}) where {C} = Term{C}(x)

function Base.convert(::Type{Term{C,T}}, x::Variable{C}) where {C,T}
    return convert(Term{C,T}, convert(Monomial{C}, x))
end
Term{C}(x::Variable{C}) where {C} = Term{C}(convert(Monomial{C}, x))
Term(x::Variable{C}) where {C} = Term{C}(x)

function MP.convert_constant(::Type{Term{C,T}}, α) where {C,T}
    return Term{C}(convert(T, α))
end
Term{C}(α::T) where {C,T} = Term{C,T}(α, Monomial{C}())

Base.broadcastable(t::Term) = Ref(t)
#(::Type{TermContainer{C}}){C}(x::Variable{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(x::Monomial{C}) = Term(x)
#(::Type{TermContainer{C}}){C}(t::TermContainer{C}) = t
#TermContainer(x::Variable) = Term(x)
#TermContainer(x::Monomial) = Term(x)
#TermContainer(t::TermContainer) = t

#Base.convert{C, T}(::Type{TermContainer{C, T}}, x::Union{Monomial{C},Variable{C}}) = Term{C, T}(x)
#Base.convert{C, T}(::Type{TermContainer{C, T}}, α::T) = Term{C, T}(α, Monomial{C}())
#Base.convert{C, S, T}(::Type{TermContainer{C, T}}, α::S) = TermContainer{C, T}(T(α))
#(::Type{TermContainer{C}}){C, T}(α::T) = TermContainer{C, T}(α)

#Base.convert{C, T}(::Type{TermContainer{C, T}}, t::Term{C}) = Term{C, T}(t)

function MA.mutable_copy(t::T) where {T<:Term}
    return T(MA.copy_if_mutable(t.α), MA.mutable_copy(t.x))
end
Base.copy(t::Term) = MA.mutable_copy(t)

MP.coefficient(t::Term) = t.α
MP.monomial(t::Term) = t.x
_vars(t) = _vars(t.x)

function MA.operate_to!(
    t::Term,
    ::typeof(*),
    t1::MP.AbstractTermLike,
    t2::MP.AbstractTermLike,
)
    MA.operate_to!(t.α, *, coefficient(t1), coefficient(t2))
    MA.operate_to!(t.x, *, monomial(t1), monomial(t2))
    return t
end
function MA.operate!(::typeof(*), t1::Term, t2::MP.AbstractTermLike)
    MA.operate!(*, t1.α, coefficient(t2))
    MA.operate!(*, t1.x, monomial(t2))
    return t1
end
function MA.operate!(::typeof(one), t::Term)
    MA.operate!(one, t.α)
    MA.operate!(zero, t.x.z)
    return t
end
