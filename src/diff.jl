function MP.differentiate(m::Monomial{C}, x::PolyVar{C}) where C
    i = findfirst(isequal(x), _vars(m))
    if (i == nothing || i == 0) || m.z[i] == 0
        zeroterm(m)
    else
        z = copy(m.z)
        z[i] -= 1
        m.z[i] * Monomial(_vars(m), z)
    end
end

function MP.differentiate(p::Polynomial{C, T}, x::PolyVar{C}) where {C, T}
    # grlex order preserved
    i = findfirst(isequal(x), _vars(p))
    S = Base.promote_op(*, T, Int)
    if i === nothing || i == 0
        zero(Polynomial{C, S})
    else
        keep = findall(z -> z[i] > 0, p.x.Z)
        Z = copy.(p.x.Z[keep])
        a = Vector{S}(uninitialized, length(keep))
        for j in 1:length(Z)
            a[j] = p.a[keep[j]] * Z[j][i]
            Z[j][i] -= 1
        end
        Polynomial(a, MonomialVector(_vars(p), Z))
    end
end
