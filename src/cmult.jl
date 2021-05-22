# Product between PolyVar and Monomial -> Monomial
function Base.:(*)(x::PolyVar{true}, y::PolyVar{true})
    if x === y
        Monomial{true}([x], [2])
    else
        Monomial{true}(x > y ? [x,y] : [y,x], [1,1])
    end
end
function multiplyvar(v::Vector{PolyVar{true}}, x::PolyVar{true}, z)
    i = findfirst(Base.Fix2(<=, x), v)
    if (i != nothing && i > 0) && v[i] == x
        copy(v), multiplyexistingvar(i, z)
    else
        j = (i == nothing || i == 0) ? length(v)+1 : i
        insertvar(v, x, j), insertvar(z, x, j)
    end
end
function Base.:(*)(x::PolyVar{true}, y::Monomial{true})
    w, z = multiplyvar(y.vars, x, y.z)
    Monomial{true}(w, z)
end
Base.:(*)(y::MonomialVector{true}, x::PolyVar{true}) = x * y
function Base.:(*)(x::PolyVar{true}, y::MonomialVector{true})
    w, Z = multiplyvar(y.vars, x, y.Z)
    MonomialVector{true}(w, Z)
end
function multdivmono!(output_variables::Vector{PolyVar{true}},
                      v::Vector{PolyVar{true}}, x::Monomial{true}, op)
    if v == x.vars
        if output_variables == v
            updatez! = (output, z) -> begin
                @. output = op(z, x.z)
                return
            end
        else
            resize!(output_variables, length(v))
            copyto!(output_variables, v)
            updatez! = (output, z) -> begin
                resize!(output, length(output_variables))
                @. output = op(z, x.z)
                return
            end
        end
    else
        maps = mergevars_to!(output_variables, [v, x.vars])
        n = length(v)
        updatez! = (output, z) -> begin
            resize!(output, length(output_variables))
            I = maps[1]; i = 1; lI = length(I)
            J = maps[2]; j = 1; lJ = length(J)
            while i <= lI || j <= lJ
                if i > lI || (j <= lJ && J[j] < I[i])
                    output[J[j]] = op(0, x.z[j])
                    j += 1
                elseif j > lJ || (i <= lI && I[i] < J[j])
                    output[I[i]] = op(z[I[i]], 0)
                    i += 1
                else
                    @assert I[i] == J[j]
                    output[I[i]] = op(z[I[i]], x.z[j])
                    i += 1
                    j += 1
                end
            end
        end
    end
    return updatez!
end
function _operate_exponents(op::F, z1::Vector{Int}, z2::Vector{Int}) where {F<:Function}
    return op.(z1, z2)
end
function _operate_exponents(op::F, Z::Vector{Vector{Int}}, z2::Vector{Int}) where {F<:Function}
    return Vector{Int}[_operate_exponents(op, z, z2) for z in Z]
end
function _operate_exponents(op::F, z1::Vector{Int}, z2::Vector{Int}, n, maps) where {F<:Function}
    newz = zeros(Int, n)
    I = maps[1]; i = 1; lI = length(I)
    J = maps[2]; j = 1; lJ = length(J)
    while i <= lI || j <= lJ
        if i > lI || (j <= lJ && J[j] < I[i])
            newz[J[j]] = op(0, z2[j])
            j += 1
        elseif j > lJ || (i <= lI && I[i] < J[j])
            newz[I[i]] = op(z1[i], 0)
            i += 1
        else
            @assert I[i] == J[j]
            newz[I[i]] = op(z1[i], z2[j])
            i += 1
            j += 1
        end
    end
    return newz
end
function _operate_exponents(op::F, Z::Vector{Vector{Int}}, z2::Vector{Int}, n, maps) where {F<:Function}
    return Vector{Int}[_operate_exponents(op, z, z2, n, maps) for z in Z]
end
function multdivmono(v::Vector{PolyVar{true}}, x::Monomial{true}, op, z)
    if v == x.vars
        w = copy(v)
        z_new = _operate_exponents(op, z, x.z)
    else
        w, maps = mergevars([v, x.vars])
        z_new = _operate_exponents(op, z, x.z, length(w), maps)
    end
    return w, z_new
end
function MP.mapexponents_to!(output::Monomial{true}, f::Function, x::Monomial{true}, y::Monomial{true})
    updatez! = multdivmono!(output.vars, x.vars, y, f)
    updatez!(output.z, x.z)
    return x
end
function MP.mapexponents!(f::Function, x::Monomial{true}, y::Monomial{true})
    updatez! = multdivmono!(x.vars, x.vars, y, f)
    updatez!(x.z, x.z)
    return x
end
function MP.mapexponents(f::Function, x::Monomial{true}, y::Monomial{true})
    w, z = multdivmono(x.vars, y, f, x.z)
    return Monomial{true}(w, z)
end
function Base.:(*)(x::Monomial{true}, y::MonomialVector{true})
    w, Z = multdivmono(y.vars, x, +, y.Z)
    return MonomialVector{true}(w, Z)
end
Base.:(*)(y::MonomialVector{true}, x::Monomial{true}) = x * y
Base.:(*)(x::Monomial{true}, y::PolyVar{true}) = y * x
