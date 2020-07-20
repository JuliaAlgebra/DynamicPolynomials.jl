MP.promote_rule_constant(::Type{T}, ::Type{<:DMonomialLike{C}}) where {C, T} = Term{C, promote_type(T, Int)}
MP.promote_rule_constant(::Type{S}, ::Type{Term{C, T}}) where {S, C, T} = Term{C, promote_type(S, T)}
MP.promote_rule_constant(::Type{S}, ::Type{<:TermPoly{C, T}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}
