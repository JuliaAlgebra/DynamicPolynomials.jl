function Base.:(*)(x::Variable{V,M}, y::Variable{V,M}) where {V<:NonCommutative,M}
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

function Base.:(*)(x::Monomial{V,M}, y::Monomial{V,M}) where {V<:NonCommutative,M}
    i = findlast(z -> z > 0, x.z)
    if i === nothing || i == 0
        return y
    end
    j = findfirst(z -> z > 0, y.z)
    if j === nothing || j == 0
        return x
    end
    if x.vars[i] == y.vars[j]
        w = [x.vars[1:i]; y.vars[j+1:end]]
        z = [x.z[1:i-1]; x.z[i] + y.z[j]; y.z[j+1:end]]
    else
        w = [x.vars[1:i]; y.vars[j:end]]
        z = [x.z[1:i]; y.z[j:end]]
    end
    return Monomial(w, z)
end

function Base.:(*)(y::MonomialVector{V,M}, x::DMonomialLike{V,M}) where {V<:NonCommutative,M}
    return MonomialVector{V,M}([yi * x for yi in y])
end
function Base.:(*)(x::DMonomialLike{V,M}, y::MonomialVector{V,M}) where {V<:NonCommutative,M}
    # The order may change
    # Example: y * [x^2, y^2] == [y^3, yx^2]
    return MonomialVector{V,M}([x * yi for yi in y])
end
