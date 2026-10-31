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

function Base.cmp(
    x::Variable{<:AnyCommutative{NameOrder}},
    y::Variable{<:AnyCommutative{NameOrder}},
)
    return cmp((name(y), y.kind), (name(x), x.kind))
end

function Base.hash(x::Variable{<:AnyCommutative{NameOrder}}, h::UInt)
    return hash((name(x), x.kind), h)
end

Base.:(==)(x::Variable, y::Variable) = iszero(cmp(x, y))
