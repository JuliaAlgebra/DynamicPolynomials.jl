function MP.promote_variables(m1::Monomial, m2::Monomial)
    if MP.variables(m1) == MP.variables(m2)
        return m1, m2
    end
    allvars, maps = mergevars([MP.variables(m1), MP.variables(m2)])
    z1 = zeros(Int, length(allvars))
    z1[maps[1]] = m1.z
    z2 = zeros(Int, length(allvars))
    z2[maps[2]] = m2.z
    return Monomial(allvars, z1), Monomial(allvars, z2)
end

# StarAlgebras promote_with_map implementations
# These are called by SA.maybe_promote(p, all_vars, map) when map is an ExponentMap,
# as part of SA.promote_bases_with_maps (defined in MultivariatePolynomials.jl).
# Guarded by isdefined since ExponentMap is not yet in released MultivariatePolynomials.

function SA.promote_with_map(
    v::Variable{V,M},
    all_vars::Vector{Variable{V,M}},
    map::MP.ExponentMap,
) where {V,M}
    new_z = map([1])
    return Monomial{V,M}(copy(all_vars), new_z), map
end

function SA.promote_with_map(
    m::Monomial{V,M},
    all_vars::Vector{Variable{V,M}},
    map::MP.ExponentMap,
) where {V,M}
    return Monomial{V,M}(copy(all_vars), map(m.z)), map
end

function SA.promote_with_map(
    p::Polynomial{V,M,T},
    all_vars::Vector{Variable{V,M}},
    map::MP.ExponentMap,
) where {V,M,T}
    new_Z = [map(z) for z in p.x.Z]
    new_x = MonomialVector{V,M}(copy(all_vars), new_Z)
    return Polynomial{V,M,T}(copy(p.a), new_x), map
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
