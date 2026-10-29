struct _NoncommutativePowers{M}
    monomial::M
end

Base.IteratorSize(::Type{<:_NoncommutativePowers}) = Base.SizeUnknown()

# Ignore padding and combine powers that become adjacent after removing zeros.
# Inlining keeps the returned variable/exponent tuples on the stack.
@inline function Base.iterate(powers::_NoncommutativePowers, i::Int = 1)
    vars, exps = MP.variables(powers.monomial), MP.exponents(powers.monomial)
    n = length(exps)
    while i <= n && iszero(exps[i])
        i += 1
    end
    i > n && return nothing
    var, exp = vars[i], exps[i]
    i += 1
    while i <= n
        if iszero(exps[i])
            i += 1
        elseif vars[i] == var
            exp += exps[i]
            i += 1
        else
            break
        end
    end
    return (var, exp), i
end

function _compare_noncommutative_powers(a, b)
    ap, bp = _NoncommutativePowers(a), _NoncommutativePowers(b)
    an, bn = iterate(ap), iterate(bp)
    acycle = bcycle = 0
    while an !== nothing && bn !== nothing
        (av, ae), ai = an
        (bv, be), bi = bn
        # The padded representation repeats the variables in decreasing order.
        # Compare its nonzero positions without constructing the zero entries.
        c = cmp(bcycle, acycle)
        iszero(c) || return c
        c = cmp(av, bv)
        iszero(c) || return c
        c = cmp(ae, be)
        iszero(c) || return c
        an, bn = iterate(ap, ai), iterate(bp, bi)
        if an !== nothing && first(first(an)) > av
            acycle += 1
        end
        if bn !== nothing && first(first(bn)) > bv
            bcycle += 1
        end
    end
    return an === nothing ? (bn === nothing ? 0 : -1) : 1
end

const _NCMonomial = DPMonomial{<:NonCommutative}

function Base.:(==)(a::_NCMonomial, b::_NCMonomial)
    return iszero(_compare_noncommutative_powers(a, b))
end
function Base.isequal(a::_NCMonomial, b::_NCMonomial)
    return iszero(_compare_noncommutative_powers(a, b))
end

function Base.cmp(a::_NCMonomial, b::_NCMonomial)
    order = MP.ordering(a)
    if order isa MP.Graded{MP.LexOrder}
        c = cmp(MP.degree(a), MP.degree(b))
        iszero(c) || return c
    elseif !(order isa MP.LexOrder)
        a, b = MP.promote_variables(a, b)
        return cmp(order, MP.exponents(a), MP.exponents(b))
    end
    return _compare_noncommutative_powers(a, b)
end

function Base.hash(m::_NCMonomial, u::UInt)
    powers = _NoncommutativePowers(m)
    next = iterate(powers)
    next === nothing && return hash(1, u)
    (var, exp), state = next
    next = iterate(powers, state)
    # The hash of a degree-one monomial must match the variable's hash.
    next === nothing && isone(exp) && return hash(var, u)
    u = hash((var, exp), u)
    while next !== nothing
        power, state = next
        u = hash(power, u)
        next = iterate(powers, state)
    end
    return u
end

# Noncommutative monomials retain the ordered variable/exponent representation:
# x^2 * y * x has variables [x, y, x] and exponents [2, 1, 1].
# Algebra elements use these complete monomials as coefficient keys.
function MP.promote_variables_with_maps(
    a::Vector{Variable{V,M}},
    b::Vector{Variable{V,M}},
) where {V<:NonCommutative,M}
    if a == b
        return (a, nothing), (b, nothing)
    end
    vars, maps = mergevars([a, b])
    ma = a == vars ? nothing : MP.ExponentMap(maps[1], length(vars))
    mb = b == vars ? nothing : MP.ExponentMap(maps[2], length(vars))
    return (vars, ma), (vars, mb)
end

function Base.:*(
    x::DPMonomial{V,M},
    y::DPMonomial{V,M},
) where {V<:NonCommutative,M}
    xv, yv = MP.variables(x), MP.variables(y)
    xe, ye = MP.exponents(x), MP.exponents(y)
    xlast = findlast(!iszero, xe)
    xlast === nothing && return y
    yfirst = findfirst(!iszero, ye)
    yfirst === nothing && return x
    ylast = findlast(!iszero, ye)
    if xv[xlast] == yv[yfirst]
        vars = [xv[1:xlast]; yv[(yfirst+1):ylast]]
        exps = [xe[1:(xlast-1)]; xe[xlast] + ye[yfirst]; ye[(yfirst+1):ylast]]
    else
        vars = [xv[1:xlast]; yv[yfirst:ylast]]
        exps = [xe[1:xlast]; ye[yfirst:ylast]]
    end
    return MP.monomial(vars, exps)
end

function Base.:^(m::DPMonomial{<:NonCommutative}, n::Integer)
    return Base.power_by_squaring(m, n)
end
