# Noncommutative monomials retain the ordered variable/exponent representation:
# x^2 * y * x has variables [x, y, x] and exponents [2, 1, 1].
# Algebra elements still use a fixed ordered sequence; products requiring a
# longer sequence are unsupported.
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
        vars = [xv[1:xlast]; yv[yfirst+1:ylast]]
        exps = [xe[1:xlast-1]; xe[xlast] + ye[yfirst]; ye[yfirst+1:ylast]]
    else
        vars = [xv[1:xlast]; yv[yfirst:ylast]]
        exps = [xe[1:xlast]; ye[yfirst:ylast]]
    end
    return MP.monomial(vars, exps)
end

function Base.:^(m::DPMonomial{<:NonCommutative}, n::Integer)
    return Base.power_by_squaring(m, n)
end
