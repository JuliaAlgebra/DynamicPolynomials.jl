# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(Variable) and one(Monomial) do not return
# a Variable and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
Base.:(^)(x::Variable{V,M}, i::Int) where {V,M} = Monomial{V,M}([x], [i])
Base.:(^)(x::Monomial{<:Commutative}, i::Int) = Monomial(x.vars, i * x.z)

myminivect(x::T, y::T) where {T} = [x, y]
function myminivect(x::S, y::T) where {S,T}
    U = promote_type(S, T)
    return [U(x), U(y)]
end

Base.:(+)(x::DMonomialLike, y::DMonomialLike) = MP.term(x) + MP.term(y)
Base.:(-)(x::DMonomialLike, y::DMonomialLike) = MP.term(x) - MP.term(y)

_getindex(p::Polynomial, i::Int) = p[i]
_getindex(t::_Term, ::Int) = t
function _plusorminus_to!(
    a::Vector{U},
    Z::Vector{Vector{Int}},
    op::Function,
    p::TermPoly{V,M},
    q::TermPoly{V,M},
    maps,
    nvars,
) where {V,M,U}
    i = j = 1
    while i <= nterms(p) || j <= nterms(q)
        z = zeros(Int, nvars)
        if j > nterms(q) || (
            i <= nterms(p) &&
            MP.monomial(_getindex(p, i)) < MP.monomial(_getindex(q, j))
        )
            t = _getindex(p, i)
            z[maps[1]] = MP.monomial(t).z
            α = MA.scaling_convert(U, MA.copy_if_mutable(MP.coefficient(t)))
            i += 1
        elseif i > nterms(p) ||
               MP.monomial(_getindex(q, j)) < MP.monomial(_getindex(p, i))
            t = _getindex(q, j)
            z[maps[2]] = MP.monomial(t).z
            α = MA.scaling_convert(U, MA.operate(op, MP.coefficient(t)))
            j += 1
        else
            t = _getindex(p, i)
            z[maps[1]] = MP.monomial(t).z
            s = _getindex(q, j)
            α = op(MP.coefficient(t), MP.coefficient(s))
            i += 1
            j += 1
        end
        push!(a, α)
        push!(Z, z)
    end
end
function plusorminus(
    p::TermPoly{V,M,S},
    q::TermPoly{V,M,T},
    op::Function,
) where {V,M,S,T}
    varsvec = [MP.variables(p), MP.variables(q)]
    allvars, maps = mergevars(varsvec)
    U = MA.promote_operation(op, S, T)
    a = U[]
    Z = Vector{Int}[]
    _plusorminus_to!(a, Z, op, p, q, maps, length(allvars))
    return Polynomial(a, MonomialVector{V,M}(allvars, Z))
end

function MA.operate!(::typeof(*), p::Polynomial, t::_Term)
    # In case `coefficient(t)` is a polynomial (e.g. in `gcd` algorithm)
    # We cannot do `MA.operate!(*, ...)`, we need `MP.right_constant_mult`
    MA.operate!(MP.right_constant_mult, p, coefficient(t))
    MA.operate!(*, p, monomial(t))
    return p
end

function MA.operate_to!(
    output::Polynomial{V,M},
    op::Union{typeof(+),typeof(-)},
    p::TermPoly{V,M},
    q::TermPoly{V,M},
) where {V,M}
    if output === p || output === q
        # Otherwise, `_plusorminus_to!` never finishes
        error(
            "Cannot call `operate_to!(output, $op, p, q)` with `output` equal to `p` or `q`, call `operate!` instead.",
        )
    end
    varsvec = [MP.variables(p), MP.variables(q)]
    allvars, maps = mergevars(varsvec)
    Future.copy!(output.x.vars, allvars)
    empty!(output.a)
    empty!(output.x.Z)
    _plusorminus_to!(output.a, output.x.Z, op, p, q, maps, length(allvars))
    return output
end

function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::Polynomial,
    q::Union{Variable,Monomial,_Term},
)
    return MA.operate!(op, p, polynomial(q))
end
const _NoVarTerm{T} = Tuple{T,Vector{Int}}
function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::Polynomial{<:NonCommutative},
    q::Polynomial{<:NonCommutative},
)
    return MA.operate_to!(p, op, MA.copy(p), q)
end

function _exponents_compare(q::Polynomial{V,M}, j, e) where {V,M}
    return MP.compare(q.x.Z[j], e, M)
end

# TODO need to check that this also works for non-commutative
function MA.operate!(
    op::Union{typeof(+),typeof(-)},
    p::Polynomial{V},
    q::Polynomial{V},
) where {V<:Commutative}
    if MP.variables(p) != MP.variables(q)
        varsvec = [MP.variables(p), MP.variables(q)]
        allvars, maps = mergevars(varsvec)
        if length(allvars) != length(MP.variables(p))
            _add_variables!(p.x, allvars, maps[1])
        end
        if length(allvars) == length(MP.variables(q))
            rhs = q
        else
            # We could avoid promoting `q` to the same variables
            # like in `plusorminus` to avoid extra allocation but it then
            # gives slower comparison. There is a tradeoff and the approach used here
            # should be better of `q` has less terms and then the same term is compared
            # many times.
            rhs = Polynomial(q.a, copy(q.x))
            _add_variables!(rhs.x, allvars, maps[2])
        end
        return MA.operate!(op, p, rhs)
    end
    get1(i) = (p.a[i], p.x.Z[i])
    function get2(i)
        return (
            MA.scaling_convert(eltype(p.a), MA.operate(op, q.a[i])),
            copy(q.x.Z[i]),
        )
    end
    function set(i, t::_NoVarTerm)
        p.a[i] = t[1]
        return p.x.Z[i] = t[2]
    end
    function push(t::_NoVarTerm)
        push!(p.a, t[1])
        return push!(p.x.Z, t[2])
    end
    compare_monomials(t::_NoVarTerm, j::Int) = _exponents_compare(q, j, t[2])
    compare_monomials(i::Int, j::Int) = compare_monomials(get1(i), j)
    combine(i::Int, j::Int) = p.a[i] = MA.operate!!(op, p.a[i], q.a[j])
    combine(t::_NoVarTerm, j::Int) = (MA.operate!!(op, t[1], q.a[j]), t[2])
    function resize(n)
        resize!(p.a, n)
        return resize!(p.x.Z, n)
    end
    # We can modify the coefficient since it's the result of `combine`.
    keep(t::_NoVarTerm) = !MA.iszero!!(t[1])
    keep(i::Int) = !MA.iszero!!(p.a[i])
    MP.polynomial_merge!(
        nterms(p),
        nterms(q),
        get1,
        get2,
        set,
        push,
        compare_monomials,
        combine,
        keep,
        resize,
    )
    return p
end

Base.:(+)(x::TermPoly{V,M}, y::TermPoly{V,M}) where {V,M} = plusorminus(x, y, +)
Base.:(-)(x::TermPoly{V,M}, y::TermPoly{V,M}) where {V,M} = plusorminus(x, y, -)
function Base.:(+)(x::TermPoly{V,M}, y::Union{Monomial,Variable}) where {V,M}
    return x + MP.term(y)
end
function Base.:(+)(x::Union{Monomial,Variable}, y::TermPoly{V,M}) where {V,M}
    return MP.term(x) + y
end

function Base.:(-)(x::TermPoly{V,M,T}, y::DMonomialLike{V,M}) where {V,M,T}
    return x - convert(MP.term_type(y, T), y)
end
function Base.:(-)(x::DMonomialLike{V,M}, y::TermPoly{V,M,T}) where {V,M,T}
    return convert(MP.term_type(x, T), x) - y
end

# `MA.operate(-, p)` redirects to `-p` as it assumes that `-p` can be modified
# through the MA API without modifying `p`. We should either copy `p.x` here
# or implement a `MA.operate(-, p)` that copies it. We choose the first option.
Base.:(-)(p::Polynomial) = Polynomial(-p.a, copy(p.x))
# TODO use that with MA v0.2
#Base.:(-)(p::Polynomial) = Polynomial(map(α -> MA.operate(-, α), p.a), copy(p.x))

include("mult.jl")
