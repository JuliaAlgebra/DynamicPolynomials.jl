function MP.antidifferentiate(m::Monomial{V,M}, x::Variable{V,M}) where {V,M}
    z = copy(m.z)
    i = findfirst(isequal(x), MP.variables(m))
    if (i === nothing || i == 0) || m.z[i] == 0
        Monomial(MP.variables(m), z) * x
    else
        z[i] += 1
        Monomial(MP.variables(m), z) / (m.z[i] + 1)
    end
end
