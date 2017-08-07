export RationalPoly
import Base.+, Base.-, Base.*, Base./

iscomm(r::Type{RationalPoly{C, S, T}}) where {C, S, T} = C

Base.convert(::Type{RationalPoly{C, S, T}}, q::RationalPoly{C, S, T}) where {C, S, T} = q
Base.convert(::Type{RationalPoly{C, S, T}}, q::RationalPoly{C, U, V}) where {C, S, T, U, V} = TermContainer{C, S}(q.num) / TermContainer{C, T}(q.den)
function Base.convert(::Type{RationalPoly{C, S, T}}, p::TermContainer{C, S}) where {C, S, T}
    p / one(TermContainer{C, T})
end
function Base.convert(::Type{RationalPoly{C, S, T}}, p::TermContainer) where {C, S, T}
    convert(RationalPoly{C, S, T}, TermContainer{C, S}(p))
end
function Base.convert(::Type{RationalPoly{C, S, T}}, p) where {C, S, T}
    Base.convert(RationalPoly{C, S, T}, TermContainer{C, S}(p))
end

(/)(r::RationalPoly, p::TermContainer) = r.num / (r.den * p)
function (/)(num::TermContainer{C, S}, den::TermContainer{C, T}) where {C, S, T}
    RationalPoly{C, S, T}(num, den)
end
function (/)(num, den::PolyType{C}) where {C}
    TermContainer{C}(num) / den
end
(/)(num::PolyType{C}, den::PolyType{C}) where {C} = TermContainer{C}(num) / TermContainer{C}(den)

# Polynomial divided by coefficient is a polynomial not a rational polynomial
(/)(num::PolyType{C}, den) where {C} = num * (1 / den)

function (+)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den + r.den*s.num) / (r.den * s.den)
end
function (+)(p::TermContainer, r::RationalPoly)
    (p*r.den + r.num) / r.den
end
(+)(r::RationalPoly, p::Polynomial) = p + r
(+)(r::RationalPoly, t::Term) = t + r
function (-)(r::RationalPoly, s::RationalPoly)
    (r.num*s.den - r.den*s.num) / (r.den * s.den)
end
(-)(p::PolyType, s::RationalPoly) = (p * s.den - s.num) / s.den
(-)(s::RationalPoly, p::PolyType) = (s.num - p * s.den) / s.den

(*)(r::RationalPoly, s::RationalPoly) = (r.num*s.num) / (r.den*s.den)
(*)(p::TermContainer, r::RationalPoly) = p == r.den ? r.num : (p * r.num) / r.den
(*)(r::RationalPoly, p::Polynomial) = p == r.den ? r.num : (r.num * p) / r.den
(*)(r::RationalPoly, t::Term)          = t == r.den ? r.num : (r.num * t) / r.den
(*)(p::PolyType, r::RationalPoly) = TermContainer(p) * r
(*)(r::RationalPoly, p::Monomial) = r * TermContainer(p)
(*)(r::RationalPoly, p::PolyVar)  = r * TermContainer(p)
(*)(α, r::RationalPoly{C}) where {C} = TermContainer{C}(α) * r
(*)(r::RationalPoly{C}, α) where {C} = r * TermContainer{C}(α)

zero(r::RationalPoly) = zero(r.num)
zero(::Type{RationalPoly{C, S, T}}) where {C, S, T} = zero(Polynomial{C, S})
