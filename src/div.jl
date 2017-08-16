function MP.divides(m1::Monomial, m2::Monomial)
    e1 = exponents(m1)
    v1 = variables(m1)
    e2 = exponents(m2)
    v2 = variables(m2)
    i = j = 1
    while i <= length(e1) && j <= length(e2)
        if v1[i] == v2[j]
            if e1[i] > e2[j]
                return false
            end
            i += 1
            j += 1
        elseif v1[i] > v2[j]
            if !iszero(e1[i])
                return false
            end
            i += 1
        else
            j += 1
        end
    end
    i > length(e1)
end
