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

function MP.antidifferentiate(p::Polynomial{V,M,T}, x::Variable{V,M}) where {V,M,T}
    i = something(findfirst(isequal(x), MP.variables(p)), 0)
    S = typeof(zero(T) * 0)
    if iszero(i)
        x * p
    else
        Z = copy.(p.x.Z)
        a = Vector{S}(undef, length(p.a))
        for j in 1:length(Z)
            a[j] = p.a[j] / (Z[j][i] + 1)
            Z[j][i] += 1
        end
        Polynomial(a, MonomialVector(MP.variables(p), Z))
    end
end
