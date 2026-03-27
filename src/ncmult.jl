function Base.:(*)(
    x::Variable{V,M},
    y::Variable{V,M},
) where {V<:NonCommutative,M}
    if x === y
        Monomial{V,M}([x], [2])
    else
        Monomial{V,M}([x, y], [1, 1])
    end
end

function multiplyvar(
    v::Vector{Variable{V,M}},
    z::Vector{Int},
    x::Variable{V,M},
) where {V,M}
    i = length(v)
    while i > 0 && z[i] == 0
        i -= 1
    end
    if i > 0 && v[i] == x
        return copy(v), multiplyexistingvar(i, z)
    else
        #   ---->
        # \  |\  |\
        #  \ | \ | \
        #   \|  \|  \
        # If v[i] > x, we wait either for a rise (v[i] > v[i-1]) or v[i] < x
        # Otherwise, we first wait for a drop and then wait for the same thing
        ndrop = 0
        if i > 0 && v[i] > x
            droplim1 = 0
            droplim2 = 1
        else
            droplim1 = 1
            droplim2 = 2
        end
        i += 1
        while i <= length(v) && v[i] != x
            if i > 1 && v[i] > v[i-1]
                ndrop += 1
            end
            if ndrop >= droplim2 || (ndrop >= droplim1 && v[i] < x)
                break
            end
            i += 1
        end

        if i <= length(v) && v[i] == x
            return copy(v), multiplyexistingvar(i, z)
        else
            return insertvar(v, x, i), insertvar(z, x, i)
        end
    end
end
function multiplyvar(
    x::Variable{<:NonCommutative},
    v::Vector{<:Variable{<:NonCommutative}},
    z::Vector{Int},
)
    i = 1
    while i <= length(v) && z[i] == 0
        i += 1
    end
    if i <= length(v) && v[i] == x
        return copy(v), multiplyexistingvar(i, z)
    else
        #   <----
        # \  |\  |\
        #  \ | \ | \
        #   \|  \|  \
        # If z[i] < x, we wait either for a drop (v[i] < v[i+1]) or v[i] > x
        # Otherwise, we first wait for a drop and then wait for the same thing
        ndrop = 0
        if i <= length(v) && v[i] < x
            droplim1 = 0
            droplim2 = 1
        else
            droplim1 = 1
            droplim2 = 2
        end
        i -= 1
        while i > 0 && v[i] != x
            if i < length(v) && v[i] < v[i+1]
                ndrop += 1
            end
            if ndrop >= droplim2 || (ndrop >= droplim1 && v[i] > x)
                break
            end
            i -= 1
        end
        if i > 0 && v[i] == x
            return copy(v), multiplyexistingvar(i, z)
        else
            return insertvar(v, x, i + 1), insertvar(z, x, i + 1)
        end
    end
end
function Base.:(*)(x::Variable{<:NonCommutative}, y::Monomial{<:NonCommutative})
    w, z = multiplyvar(x, y.vars, y.z)
    return Monomial(w, z)
end
function Base.:(*)(y::Monomial{<:NonCommutative}, x::Variable{<:NonCommutative})
    w, z = multiplyvar(y.vars, y.z, x)
    return Monomial(w, z)
end

function Base.:(*)(
    x::Monomial{V,M},
    y::Monomial{V,M},
) where {V<:NonCommutative,M}
    xj = findlast(!iszero, x.z)
    if xj === nothing || xj == 0
        return y
    end
    yi = findfirst(!iszero, y.z)
    if yi === nothing || yi == 0
        return x
    end
    yj = findlast(!iszero, y.z)
    if x.vars[xj] == y.vars[yi]
        w = [x.vars[1:xj]; y.vars[yi+1:yj]]
        z = [x.z[1:xj-1]; x.z[xj] + y.z[yi]; y.z[yi+1:yj]]
    else
        w = [x.vars[1:xj]; y.vars[yi:yj]]
        z = [x.z[1:xj]; y.z[yi:yj]]
    end
    return Monomial(w, z)
end

function Base.:(*)(
    y::MonomialVector{V,M},
    x::DMonomialLike{V,M},
) where {V<:NonCommutative,M}
    return MonomialVector{V,M}([yi * x for yi in y])
end
function Base.:(*)(
    x::DMonomialLike{V,M},
    y::MonomialVector{V,M},
) where {V<:NonCommutative,M}
    # The order may change
    # Example: y * [x^2, y^2] == [y^3, yx^2]
    return MonomialVector{V,M}([x * yi for yi in y])
end
