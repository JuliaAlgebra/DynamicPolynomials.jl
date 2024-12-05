function MP.promote_rule_constant(
    ::Type{Any},
    ::Type{<:DMonomialLike{V,M}},
) where {V,M}
    return Any
end
function MP.promote_rule_constant(
    ::Type{T},
    ::Type{<:DMonomialLike{V,M}},
) where {V,M,T}
    return _Term{V,M,promote_type(T, Int)}
end
function MP.promote_rule_constant(
    ::Type{S},
    ::Type{_Term{V,M,T}},
) where {S,V,M,T}
    return _Term{V,M,promote_type(S, T)}
end
function MP.promote_rule_constant(
    ::Type{S},
    ::Type{<:TermPoly{V,M,T}},
) where {S,V,M,T}
    return Polynomial{V,M,promote_type(S, T)}
end
MP.promote_rule_constant(::Type, ::Type{_Term{V,M}}) where {V,M} = Any
MP.promote_rule_constant(::Type, ::Type{Polynomial{V,M}}) where {V,M} = Any
function Base.promote_rule(
    ::Type{_Term{V,M}},
    ::Type{_Term{V,M,T}},
) where {V,M,T}
    return _Term{V,M}
end
function Base.promote_rule(
    ::Type{_Term{V,M,T}},
    ::Type{_Term{V,M}},
) where {V,M,T}
    return _Term{V,M}
end
function Base.promote_rule(
    ::Type{_Term{V,M}},
    ::Type{<:DMonomialLike{V,M}},
) where {V,M}
    return _Term{V,M}
end
function Base.promote_rule(
    ::Type{<:DMonomialLike{V,M}},
    ::Type{_Term{V,M}},
) where {V,M}
    return _Term{V,M}
end

function Base.convert(::Type{_Term{V,M}}, m::DMonomialLike{V,M}) where {V,M}
    return convert(_Term{V,M,Int}, m)
end
function Base.convert(
    ::Type{Polynomial{V,M}},
    t::Union{TermPoly{V,M},DMonomialLike{V,M}},
) where {V,M}
    return MP.polynomial(t)
end
