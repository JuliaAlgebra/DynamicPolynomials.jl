import Base.promote_rule

# Promotion with PolyVar and Monomial
promote_rule{C}(::Type{Monomial{C}}, ::Type{PolyVar{C}}) = Monomial{C}
promote_rule{C}(::Type{PolyVar{C}}, ::Type{Monomial{C}}) = Monomial{C}
#promote_rule{S<:Union{Monomial, PolyVar}, T<:Union{Monomial, PolyVar}}(::Type{S}, ::Type{T}) = Monomial{iscomm{S}}

promote_rule{S, C, T}(::Type{<:AbstractPolynomialLike{S}}, ::Type{<:TermPoly{C, T}}) = Polynomial{C, promote_type(S, T)}
promote_rule{S, C, T}(::Type{<:TermPoly{C, T}}, ::Type{<:AbstractPolynomialLike{S}}) = Polynomial{C, promote_type(S, T)}

promote_rule{S, C, T}(::Type{S}, ::Type{<:TermPoly{C, T}}) = Polynomial{C, promote_type(S, T)}
promote_rule{S, C, T}(::Type{<:TermPoly{C, T}}, ::Type{S}) = Polynomial{C, promote_type(S, T)}
promote_rule{C, T}(::Type{T}, ::Type{<:DMonomialLike{C}}) = Term{C, promote_type(T, Int)}
promote_rule{C, T}(::Type{<:DMonomialLike{C}}, ::Type{T}) = Term{C, promote_type(T, Int)}

# Promotion with Term
promote_rule{C,S,T}(::Type{Term{C, S}}, ::Type{Term{C, T}}) = Term{C, promote_type(S, T)}
promote_rule{S,C,T}(::Type{S}, ::Type{Term{C, T}}) = Term{C, promote_type(S, T)}
promote_rule{S,C,T}(::Type{Term{C, T}}, ::Type{S}) = Term{C, promote_type(S, T)}
promote_rule{C,T}(::Type{Monomial{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{Monomial{C}}) = Term{C, T}
promote_rule{C,T}(::Type{PolyVar{C}}, ::Type{Term{C, T}}) = Term{C, T}
promote_rule{C,T}(::Type{Term{C, T}}, ::Type{PolyVar{C}}) = Term{C, T}

# Promotion with Polynomial
promote_rule{C, S, T}(::Type{Polynomial{C, S}}, ::Type{Term{C, T}}) = Polynomial{C, promote_type(S, T)}
promote_rule{C, S, T}(::Type{Term{C, T}}, ::Type{Polynomial{C, S}}) = Polynomial{C, promote_type(S, T)}

promote_rule{C, S, T}(::Type{<:TermPoly{C, S}}, ::Type{<:TermPoly{C, T}}) = Polynomial{C, promote_type(S, T)}
promote_rule{C, T}(::Type{<:TermPoly{C, T}}, ::Type{<:DMonomialLike{C}}) = Polynomial{C, promote_type(Int, T)}
promote_rule{C, T}(::Type{<:DMonomialLike{C}}, ::Type{<:TermPoly{C, T}}) = Polynomial{C, promote_type(Int, T)}
