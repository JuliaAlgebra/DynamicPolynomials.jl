function MP.differentiate(m::Monomial{C}, x::PolyVar{C}) where C
    i = findfirst(_vars(m), x)
    if i == 0 || m.z[i] == 0
        zeroterm(m)
    else
        z = copy(m.z)
        z[i] -= 1
        m.z[i] * Monomial(_vars(m), z)
    end
end

function MP.differentiate(p::Polynomial{C, T}, x::PolyVar{C}) where {C, T}
    # grlex order preserved
    i = findfirst(_vars(p), x)
    S = Base.promote_op(*, T, Int)
    if i == 0
        zero(Polynomial{C, S})
    else
        keep = find([z[i] > 0 for z in p.x.Z])
        Z = [copy(p.x.Z[i]) for i in keep]
        a = Vector{S}(length(keep))
        for j in 1:length(Z)
            a[j] = p.a[keep[j]] * Z[j][i]
            Z[j][i] -= 1
        end
        Polynomial(a, MonomialVector(_vars(p), Z))
    end
end
