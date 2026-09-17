import Base.==

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
