import Base.==

# Comparison of Variable vectors
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

Base.:(==)(x::Variable, y::Variable) = iszero(cmp(x, y))
