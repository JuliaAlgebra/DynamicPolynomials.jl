# Product between PolyVar and Monomial -> Monomial
function (*)(x::PolyVar{true}, y::PolyVar{true})
    if x === y
        Monomial{true}([x], [2])
    else
        Monomial{true}(x > y ? [x,y] : [y,x], [1,1])
    end
end
function multiplyvar(v::Vector{PolyVar{true}}, x::PolyVar{true})
    i = findfirst(w->w <= x, v)
    if i > 0 && v[i] == x
        multiplyexistingvar(v, x, i)
    else
        insertvar(v, x, i == 0 ? length(v)+1 : i)
    end
end
function (*)(x::PolyVar{true}, y::Monomial{true})
    w, updatez = multiplyvar(y.vars, x)
    Monomial{true}(w, updatez(y.z))
end
(*)(y::MonomialVector{true}, x::PolyVar{true}) = x * y
function (*)(x::PolyVar{true}, y::MonomialVector{true})
    w, updatez = multiplyvar(y.vars, x)
    MonomialVector{true}(w, updatez.(y.Z))
end
function multdivmono(v::Vector{PolyVar{true}}, x::Monomial{true}, op)
    if v == x.vars
        # /!\ no copy done here for efficiency, do not mess up with vars
        w = v
        updatez = z -> op(z, x.z)
    else
        w, maps = myunion([v, x.vars])
        updatez = z -> begin
            newz = zeros(Int, length(w))
            newz[maps[1]] += z
            newz[maps[2]] = op(newz[maps[2]], x.z)
            newz
        end
    end
    w, updatez
end
function (*)(x::Monomial{true}, y::Monomial{true})
    w, updatez = multdivmono(y.vars, x, +)
    Monomial{true}(w, updatez(y.z))
end
function (*)(x::Monomial{true}, y::MonomialVector{true})
    w, updatez = multdivmono(y.vars, x, +)
    MonomialVector{true}(w, updatez.(y.Z))
end
(*)(y::MonomialVector{true}, x::Monomial{true}) = x * y
(*)(x::Monomial{true}, y::PolyVar{true}) = y * x
