function MP.differentiate(m::Monomial{C}, x::PolyVar{C}) where C
    i = findfirst(isequal(x), _vars(m))
    if (i === nothing || i == 0) || m.z[i] == 0
        zero_term(m)
    else
        z = copy(m.z)
        z[i] -= 1
        m.z[i] * Monomial(_vars(m), z)
    end
end

function MP.differentiate(p::Polynomial{C, T}, x::PolyVar{C}) where {C, T}
    # grlex order preserved
    i = something(findfirst(isequal(x), _vars(p)), 0)
    S = typeof(zero(T) * 0)
    if iszero(i)
        zero(Polynomial{C, S})
    else
        keep = findall(z -> z[i] > 0, p.x.Z)
        Z = copy.(p.x.Z[keep])
        a = Vector{S}(undef, length(keep))
        for j in 1:length(Z)
            a[j] = p.a[keep[j]] * Z[j][i]
            Z[j][i] -= 1
        end
        Polynomial(a, MonomialVector(_vars(p), Z))
    end
end

function MP.antidifferentiate(m::Monomial{C}, x::PolyVar{C}) where {C}
    i = findfirst(isequal(x), _vars(m))
    if (i === nothing || i == 0) || m.z[i] == 0
        # Insert `x` in the monomial
        return m * x
    else
        z = copy(m.z)
        z[i] += 1
        (1 // (m.z[i] + 1)) * Monomial(_vars(m), z)
    end
end

function MP.antidifferentiate(p::Polynomial{C,T}, x::PolyVar{C}) where {C,T}
    # Increase by one the order of `x` in all monomials of `p`
    q = p * x

    # Find index of `x` in monovec of `q`
    i = something(findfirst(isequal(x), _vars(q)), 0)

    # Divide the coef of each monomial by `1 + deg(p, x)`
    # (i.e by `deg(q, x)`)
    _x = MP.monomials(q)
    for j in eachindex(_x)
        q.a[j] //= _x.Z[j][i]
    end

    return q
end