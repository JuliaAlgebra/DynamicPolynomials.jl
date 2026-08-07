function MP.differentiate(m::Monomial{V,M}, x::Variable{V,M}) where {V,M}
    i = findfirst(isequal(x), MP.variables(m))
    if (i === nothing || i == 0) || m.z[i] == 0
        zero_term(m)
    else
        z = copy(m.z)
        z[i] -= 1
        m.z[i] * Monomial(MP.variables(m), z)
    end
end
