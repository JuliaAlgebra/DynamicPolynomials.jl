# StarAlgebras promote_with_map for DP Variables
function SA.promote_with_map(
    v::Variable{V,M},
    all_vars::Vector{Variable{V,M}},
    map::MP.ExponentMap,
) where {V,M}
    new_z = map([1])
    return MP.Polynomial(MP.Variables{MP.Monomial}(copy(all_vars)), new_z), map
end
