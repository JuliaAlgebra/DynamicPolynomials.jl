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

# StarAlgebras promote_with_map implementations for Variables and Monomials
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

function MP.promote_rule_constant(
    ::Type{T},
    ::Type{<:DMonomialLike{V,M}},
) where {V,M,T}
    return _Term{V,M,promote_type(T, Int)}
end
