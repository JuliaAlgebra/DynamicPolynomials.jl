# Product between PolyVar and Monomial -> Monomial
function Base.:(*)(x::PolyVar{true}, y::PolyVar{true})
    if x === y
        Monomial{true}([x], [2])
    else
        Monomial{true}(x > y ? [x,y] : [y,x], [1,1])
    end
end
function multiplyvar(v::Vector{PolyVar{true}}, x::PolyVar{true})
    i = findfirst(w->w <= x, v)
    if (i != nothing && i > 0) && v[i] == x
        multiplyexistingvar(v, x, i)
    else
        insertvar(v, x, (i == nothing || i == 0) ? length(v)+1 : i)
    end
end
function Base.:(*)(x::PolyVar{true}, y::Monomial{true})
    w, updatez = multiplyvar(y.vars, x)
    Monomial{true}(w, updatez(y.z))
end
Base.:(*)(y::MonomialVector{true}, x::PolyVar{true}) = x * y
function Base.:(*)(x::PolyVar{true}, y::MonomialVector{true})
    w, updatez = multiplyvar(y.vars, x)
    MonomialVector{true}(w, updatez.(y.Z))
end
function multdivmono(v::Vector{PolyVar{true}}, x::Monomial{true}, op)
    if v == x.vars
        # /!\ no copy done here for efficiency, do not mess up with vars
        w = v
        updatez = z -> op.(z, x.z)
    else
        w, maps = mergevars([v, x.vars])
        updatez = z -> begin
            newz = zeros(Int, length(w))
            I = maps[1]; i = 1; lI = length(I)
            J = maps[2]; j = 1; lJ = length(J)
            while i <= lI || j <= lJ
                if i > lI || (j <= lJ && J[j] < I[i])
                    newz[J[j]] = op(0, x.z[j])
                    j += 1
                elseif j > lJ || (i <= lI && I[i] < J[j])
                    newz[I[i]] = op(z[i], 0)
                    i += 1
                else
                    @assert I[i] == J[j]
                    newz[I[i]] = op(z[i], x.z[j])
                    i += 1
                    j += 1
                end
            end
            newz
        end
    end
    w, updatez
end
function MP.mapexponents(f, x::Monomial{true}, y::Monomial{true})
    w, updatez = multdivmono(x.vars, y, f)
    Monomial{true}(w, updatez(x.z))
end
function Base.:(*)(x::Monomial{true}, y::MonomialVector{true})
    w, updatez = multdivmono(y.vars, x, +)
    MonomialVector{true}(w, updatez.(y.Z))
end
Base.:(*)(y::MonomialVector{true}, x::Monomial{true}) = x * y
Base.:(*)(x::Monomial{true}, y::PolyVar{true}) = y * x
