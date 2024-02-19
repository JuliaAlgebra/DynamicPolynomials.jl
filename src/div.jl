function MP.divides(m1::Monomial{<:NonCommutative}, m2::Monomial{<:NonCommutative})
    error("Not implemented yet")
end
function MP.divides(m1::Monomial{<:Commutative}, m2::Monomial{<:Commutative})
    e1 = exponents(m1)
    v1 = variables(m1)
    e2 = exponents(m2)
    v2 = variables(m2)
    i = 1
    lI = length(e1)
    j = 1
    lJ = length(e2)
    while i <= lI || j <= lJ
        if i > lI || (j <= lJ && v2[j] > v1[i])
            j += 1
        elseif j > lJ || (i <= lI && v1[i] > v2[j])
            iszero(e1[i]) || return false
            i += 1
        else
            @assert v1[i] == v2[j]
            e1[i] <= e2[j] || return false
            i += 1
            j += 1
        end
    end
    return true
end
