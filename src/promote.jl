MP.promote_rule_constant(::Type{T}, ::Type{<:DMonomialLike{C}}) where {C, T} = Term{C, promote_type(T, Int)}
MP.promote_rule_constant(::Type{S}, ::Type{Term{C, T}}) where {S, C, T} = Term{C, promote_type(S, T)}
MP.promote_rule_constant(::Type{S}, ::Type{<:TermPoly{C, T}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}
MP.promote_rule_constant(::Type, ::Type{Term{C}}) where {C} = Any
MP.promote_rule_constant(::Type, ::Type{Polynomial{C}}) where {C} = Any
Base.promote_rule(::Type{Term{C}}, ::Type{Term{C, T}}) where {C, T} = Term{C}
Base.promote_rule(::Type{Term{C, T}}, ::Type{Term{C}}) where {C, T} = Term{C}
Base.promote_rule(::Type{Term{C}}, ::Type{<:DMonomialLike{C}}) where {C} = Term{C}
Base.promote_rule(::Type{<:DMonomialLike{C}}, ::Type{Term{C}}) where {C} = Term{C}

Base.convert(::Type{Term{C}}, m::DMonomialLike{C}) where {C} = convert(Term{C,Int}, m)
Base.convert(::Type{Polynomial{C}}, t::Union{TermPoly{C},DMonomialLike{C}}) where {C} = MP.polynomial(t)
