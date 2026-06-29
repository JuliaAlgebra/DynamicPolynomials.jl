import Base.==

# TODO This should be in Base with T instead of Variable{V,M}.
# See https://github.com/blegat/MultivariatePolynomials.jl/issues/3
function Base.:(==)(x::Vector{Variable{V,M}}, y::Vector{Variable{V,M}}) where {V,M}
    if length(x) != length(y)
        false
    else
        #for (xi, yi) in zip(x, y)
        for i in 1:length(x)
            if x[i] != y[i]
                return false
            end
        end
        true
    end
end

# Comparison of Variable

const AnyCommutative{O} = Union{Commutative{O},NonCommutative{O}}

function Base.cmp(
    x::Variable{<:AnyCommutative{CreationOrder}},
    y::Variable{<:AnyCommutative{CreationOrder}},
)
    if x.variable_order.order.id == y.variable_order.order.id
        return cmp(y.kind, x.kind)
    else
        return cmp(y.variable_order.order.id, x.variable_order.order.id)
    end
end

# TODO remove
Base.:(==)(x::Variable, y::Variable) = iszero(cmp(x, y))
Base.:(==)(x::Monomial, y::Monomial) = iszero(cmp(x, y))

# Comparison of MonomialVector
function (==)(x::MonomialVector{V,M}, y::MonomialVector{V,M}) where {V,M}
    if length(x.Z) != length(y.Z)
        return false
    end
    allvars, maps = mergevars([MP.variables(x), MP.variables(y)])
    # Should be sorted in the same order since the non-common
    # polyvar should have exponent 0
    for (a, b) in zip(x.Z, y.Z)
        A = zeros(Int, length(allvars))
        B = zeros(Int, length(allvars))
        A[maps[1]] = a
        B[maps[2]] = b
        if A != B
            return false
        end
    end
    return true
end
(==)(mv::AbstractVector, x::MonomialVector) = monomial_vector(mv) == x
(==)(x::MonomialVector, mv::AbstractVector) = x == monomial_vector(mv)

# Comparison of Term
function _compare(p::Polynomial{V,M}, q::Polynomial{V,M}, comparator) where {V,M}
    # terms should be sorted and without zeros
    if MP.nterms(p) != MP.nterms(q)
        return false
    end
    for i in eachindex(p.a)
        if !comparator(p.x[i], q.x[i])
            # There should not be zero terms
            @assert p.a[i] != 0
            @assert q.a[i] != 0
            return false
        end
        if !comparator(p.a[i], q.a[i])
            return false
        end
    end
    return true
end

(==)(p::Polynomial{V, M}, q::Polynomial{V, M}) where {V, M} = _compare(p, q, (==))
Base.isequal(p::Polynomial{V, M}, q::Polynomial{V, M}) where {V, M} = _compare(p, q, isequal)

function Base.isapprox(
    p::Polynomial{V,M,S},
    q::Polynomial{V,M,T};
    rtol::Real = Base.rtoldefault(S, T, 0),
    atol::Real = 0,
    ztol::Real = iszero(atol) ? Base.rtoldefault(S, T, 0) : atol,
) where {V,M,S,T}
    i = j = 1
    while i <= length(p.x) || j <= length(q.x)
        if i > length(p.x) || (j <= length(q.x) && q.x[j] < p.x[i])
            if !isapproxzero(q.a[j], ztol = ztol)
                return false
            end
            j += 1
        elseif j > length(q.x) || p.x[i] < q.x[j]
            if !isapproxzero(p.a[i], ztol = ztol)
                return false
            end
            i += 1
        else
            if !isapprox(p.a[i], q.a[j], rtol = rtol, atol = atol)
                return false
            end
            i += 1
            j += 1
        end
    end
    return true
end
