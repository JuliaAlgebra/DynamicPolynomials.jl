function MP.promote_rule_constant(
    ::Type{T},
    ::Type{<:DMonomialLike{C}},
) where {C,T}
    return Term{C,promote_type(T, Int)}
end
function MP.promote_rule_constant(::Type{S}, ::Type{Term{C,T}}) where {S,C,T}
    return Term{C,promote_type(S, T)}
end
function MP.promote_rule_constant(
    ::Type{S},
    ::Type{<:TermPoly{C,T}},
) where {S,C,T}
    return Polynomial{C,promote_type(S, T)}
end
MP.promote_rule_constant(::Type, ::Type{Term{C}}) where {C} = Any
MP.promote_rule_constant(::Type, ::Type{Polynomial{C}}) where {C} = Any
Base.promote_rule(::Type{Term{C}}, ::Type{Term{C,T}}) where {C,T} = Term{C}
Base.promote_rule(::Type{Term{C,T}}, ::Type{Term{C}}) where {C,T} = Term{C}
function Base.promote_rule(
    ::Type{Term{C}},
    ::Type{<:DMonomialLike{C}},
) where {C}
    return Term{C}
end
function Base.promote_rule(
    ::Type{<:DMonomialLike{C}},
    ::Type{Term{C}},
) where {C}
    return Term{C}
end

function Base.convert(::Type{Term{C}}, m::DMonomialLike{C}) where {C}
    return convert(Term{C,Int}, m)
end
function Base.convert(
    ::Type{Polynomial{C}},
    t::Union{TermPoly{C},DMonomialLike{C}},
) where {C}
    return MP.polynomial(t)
end
