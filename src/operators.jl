# Monomial power
Base.:(^)(x::Variable{V,M}, i::Int) where {V,M} = Monomial{V,M}([x], [i])
Base.:(^)(x::Monomial{<:Commutative}, i::Int) = Monomial(copy(x.vars), i * x.z)

# Monomial arithmetic (delegates to MP term arithmetic)
Base.:(+)(x::DMonomialLike, y::DMonomialLike) = MP.term(x) + MP.term(y)
Base.:(-)(x::DMonomialLike, y::DMonomialLike) = MP.term(x) - MP.term(y)
