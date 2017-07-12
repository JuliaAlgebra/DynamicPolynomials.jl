import MultivariatePolynomials.differentiate

differentiate(p::PolyVar, x)  = differentiate(term(p), x)
differentiate(p::Monomial, x) = differentiate(term(p), x)

function differentiate{C, T}(t::Term{C, T}, x::PolyVar{C})
    i = findfirst(_vars(t), x)
    if i == 0 || t.x.z[i] == 0
        S = Base.promote_op(*, T, Int)
        zero(Term{C, S})
    else
        z = copy(t.x.z)
        z[i] -= 1
        Term(t.α * t.x.z[i], Monomial(_vars(t), z))
    end
end

function differentiate{C, T}(p::Polynomial{C, T}, x::PolyVar{C})
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
