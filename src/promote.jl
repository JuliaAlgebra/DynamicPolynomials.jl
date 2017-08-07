# Promotion with PolyVar and Monomial
Base.promote_rule(::Type{PolyVar{C}}, ::Type{PolyVar{C}}) where {C} = PolyVar{C}
Base.promote_rule(::Type{<:DMonomialLike{C}}, ::Type{<:DMonomialLike{C}}) where {C} = Monomial{C}
#promote_rule{S<:Union{Monomial, PolyVar}, T<:Union{Monomial, PolyVar}}(::Type{S}, ::Type{T}) = Monomial{iscomm{S}}

Base.promote_rule(::Type{Polynomial{C, S}}, ::Type{Polynomial{C, T}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}
Base.promote_rule(::Type{<:AbstractPolynomialLike{S}}, ::Type{Polynomial{C, T}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}
Base.promote_rule(::Type{Polynomial{C, T}}, ::Type{<:AbstractPolynomialLike{S}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}
Base.promote_rule(::Type{<:DMonomialLike{S}}, ::Type{Polynomial{C, T}}) where {S, C, T} = Polynomial{C, promote_type(T, Int)}
Base.promote_rule(::Type{Polynomial{C, T}}, ::Type{<:DMonomialLike{S}}) where {S, C, T} = Polynomial{C, promote_type(T, Int)}

MP.promote_rule_constant(::Type{T}, ::Type{<:DMonomialLike{C}}) where {C, T} = Term{C, promote_type(T, Int)}
MP.promote_rule_constant(::Type{S}, ::Type{Term{C, T}}) where {S, C, T} = Term{C, promote_type(S, T)}
MP.promote_rule_constant(::Type{S}, ::Type{<:TermPoly{C, T}}) where {S, C, T} = Polynomial{C, promote_type(S, T)}

# Promotion with Term
Base.promote_rule(::Type{Term{C, S}}, ::Type{Term{C, T}}) where {C,S,T} = Term{C, promote_type(S, T)}
Base.promote_rule(::Type{<:DMonomialLike{C}}, ::Type{Term{C, T}}) where {C,T} = Term{C, promote_type(T, Int)}
Base.promote_rule(::Type{Term{C, T}}, ::Type{<:DMonomialLike{C}}) where {C,T} = Term{C, promote_type(T, Int)}

# Promotion with Polynomial
#Base.promote_rule(::Type{Polynomial{C, S}}, ::Type{Term{C, T}}) where {C, S, T} = Polynomial{C, promote_type(S, T)}
#Base.promote_rule(::Type{Term{C, T}}, ::Type{Polynomial{C, S}}) where {C, S, T} = Polynomial{C, promote_type(S, T)}

#Base.promote_rule(::Type{<:TermPoly{C, S}}, ::Type{<:TermPoly{C, T}}) where {C, S, T} = Polynomial{C, promote_type(S, T)}
#Base.promote_rule(::Type{<:TermPoly{C, T}}, ::Type{<:DMonomialLike{C}}) where {C, T} = Polynomial{C, promote_type(T, Int)}
#Base.promote_rule(::Type{<:DMonomialLike{C}}, ::Type{<:TermPoly{C, T}}) where {C, T} = Polynomial{C, promote_type(T, Int)}
