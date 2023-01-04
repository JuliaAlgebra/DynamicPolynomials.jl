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
    if (i !== nothing && i > 0) && v[i] == x
        copy(v), multiplyexistingvar(i, z)
    else
        j = (i === nothing || i == 0) ? length(v)+1 : i
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
function _operate_exponents_to!(output::Vector{Int}, op::F, z1::Vector{Int}, z2::Vector{Int}) where {F<:Function}
    @. output = op(z1, z2)
    return
end
function _operate_exponents_to!(output::Vector{Vector{Int}}, op::F, z1::Vector{Vector{Int}}, z2::Vector{Int}) where {F<:Function}
    for i in eachindex(output)
        _operate_exponents_to!(output[i], op, z1[i], z2)
    end
    return
end
function _operate_exponents_to!(output::Vector{Int}, op::F, z1::Vector{Int}, z2::Vector{Int}, maps) where {F<:Function}
    I = maps[1]; i = 1; lI = length(I)
    J = maps[2]; j = 1; lJ = length(J)
    while i <= lI || j <= lJ
        if i > lI || (j <= lJ && J[j] < I[i])
            output[J[j]] = op(0, z2[j])
            j += 1
        elseif j > lJ || (i <= lI && I[i] < J[j])
            output[I[i]] = op(z1[i], 0)
            i += 1
        else
            @assert I[i] == J[j]
            output[I[i]] = op(z1[i], z2[j])
            i += 1
            j += 1
        end
    end
    return
end
# Not used yet
#function _operate_exponents!(op::F, Z::Vector{Vector{Int}}, z2::Vector{Int}, args::Vararg{Any,N}) where {F<:Function,N}
#    return Vector{Int}[_operate_exponents!(op, z, z2, args...) for z in Z]
#end
function _multdivmono!(output, output_variables::Vector{PolyVar{true}},
                      v::Vector{PolyVar{true}}, x::Monomial{true}, op, z)
    if v == x.vars
        if output_variables != v
            resize!(output_variables, length(v))
            copyto!(output_variables, v)
            resize!(output, length(output_variables))
        end
        _operate_exponents_to!(output, op, z, x.z)
    else
        maps = mergevars_to!(output_variables, [v, x.vars])
        resize!(output, length(output_variables))
        _operate_exponents_to!(output, op, z, x.z, maps)
    end
    return
end
function _operate_exponents(op::F, z1::Vector{Int}, z2::Vector{Int}) where {F<:Function}
    return op.(z1, z2)
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
function _operate_exponents(op::F, Z::Vector{Vector{Int}}, z2::Vector{Int}, args::Vararg{Any,N}) where {F<:Function,N}
    return Vector{Int}[_operate_exponents(op, z, z2, args...) for z in Z]
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
    if x.vars == y.vars
        if output.vars != x.vars
            n = length(x.vars)
            resize!(output.vars, n)
            copyto!(output.vars, x.vars)
            resize!(output.z, n)
        end
        _operate_exponents_to!(output.z, f, x.z, y.z)
    else
        _multdivmono!(output.z, output.vars, x.vars, y, f, x.z)
    end
    return output
end
function MP.mapexponents!(f::Function, x::Monomial{true}, y::Monomial{true})
    if x.vars == y.vars
        _operate_exponents_to!(x.z, f, x.z, y.z)
    else
        _multdivmono!(x.z, x.vars, copy(x.vars), y, f, copy(x.z))
    end
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
function MA.operate!(::typeof(*), x::MonomialVector{true}, y::Monomial{true})
    _multdivmono!(x.Z, x.vars, copy(x.vars), y, +, copy(x.Z))
    return x
end
Base.:(*)(y::MonomialVector{true}, x::Monomial{true}) = x * y
Base.:(*)(x::Monomial{true}, y::PolyVar{true}) = y * x
