import Base.==

# TODO This should be in Base with T instead of Variable{V,M}.
# See https://github.com/blegat/MultivariatePolynomials.jl/issues/3
function Base.:(==)(x::Vector{Variable{V,M}}, y::Vector{Variable{V,M}}) where {V,M}
    if length(x) != length(y)
        false
    else
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
