# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(PolyVar) and one(Monomial) do not return
# a PolyVar and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
Base.:(^)(x::PolyVar{C}, i::Int) where {C} = Monomial{C}([x], [i])
Base.:(^)(x::Monomial{true}, i::Int) = Monomial{true}(x.vars, i*x.z)

myminivect(x::T, y::T) where {T} = [x, y]
function myminivect(x::S, y::T) where {S,T}
    U = promote_type(S, T)
    [U(x), U(y)]
end

Base.:(+)(x::DMonomialLike, y::DMonomialLike) = Term(x) + Term(y)
Base.:(-)(x::DMonomialLike, y::DMonomialLike) = Term(x) - Term(y)

# `+(::UniformScaling)` is not defined...
_unary(::typeof(+), x) = MA.copy_if_mutable(x)
_unary(::typeof(-), x) = -x

_getindex(p::Polynomial, i) = p[i]
_getindex(t::Term, i) = t
function _plusorminus_to!(a::Vector{U}, Z::Vector{Vector{Int}}, op::Function, p::TermPoly{C}, q::TermPoly{C}, maps, nvars) where {C, U}
    i = j = 1
    while i <= nterms(p) || j <= nterms(q)
        z = zeros(Int, nvars)
        if j > nterms(q) || (i <= nterms(p) && _getindex(p, i).x > _getindex(q, j).x)
            t = _getindex(p, i)
            z[maps[1]] = t.x.z
            α = MA.scaling_convert(U, MA.copy_if_mutable(t.α))
            i += 1
        elseif i > nterms(p) || _getindex(q, j).x > _getindex(p, i).x
            t = _getindex(q, j)
            z[maps[2]] = t.x.z
            α = MA.scaling_convert(U, _unary(op, t.α))
            j += 1
        else
            t = _getindex(p, i)
            z[maps[1]] = t.x.z
            s = _getindex(q, j)
            α = op(t.α, s.α)
            i += 1
            j += 1
        end
        push!(a, α)
        push!(Z, z)
    end
end
function plusorminus(p::TermPoly{C, S}, q::TermPoly{C, T}, op::Function) where {C, S, T}
    varsvec = [_vars(p), _vars(q)]
    allvars, maps = mergevars(varsvec)
    U = MA.promote_operation(op, S, T)
    a = U[]
    Z = Vector{Int}[]
    _plusorminus_to!(a, Z, op, p, q, maps, length(allvars))
    Polynomial(a, MonomialVector{C}(allvars, Z))
end

function MA.mutable_operate_to!(output::Polynomial{C}, op::Union{typeof(+), typeof(-)},
                                p::TermPoly{C}, q::TermPoly{C}) where C
    if output === p || output === q
        # Otherwise, `_plusorminus_to!` never finishes
        error("Cannot call `mutable_operate_to!(output, $op, p, q)` with `output` equal to `p` or `q`, call `mutable_operate!` instead.")
    end
    varsvec = [_vars(p), _vars(q)]
    allvars, maps = mergevars(varsvec)
    Future.copy!(output.x.vars, allvars)
    empty!(output.a)
    empty!(output.x.Z)
    _plusorminus_to!(output.a, output.x.Z, op, p, q, maps, length(allvars))
    return output
end

function MA.mutable_operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial,
                             q::Union{PolyVar, Monomial, Term})
    return MA.mutable_operate!(op, p, polynomial(q))
end
const _NoVarTerm{T} = Tuple{T, Vector{Int}}
function MA.mutable_operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial{false},
                             q::Polynomial{false})
    return MA.mutable_operate_to!(p, op, MA.mutable_copy(p), q)
end
# TODO need to check that this also works for non-commutative
function MA.mutable_operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial{true},
                             q::Polynomial{true})
    if _vars(p) != _vars(q)
        varsvec = [_vars(p), _vars(q)]
        allvars, maps = mergevars(varsvec)
        _add_variables!(p.x, allvars, maps[1])
        # We could avoid promoting `q` to the same variables
        # like in `plusorminus` to avoid extra allocation but it then
        # gives slower comparison. There is a tradeoff and the approach used here
        # should be better of `q` has less terms and then the same term is compared
        # many times.
        rhs = Polynomial(q.a, copy(q.x))
        _add_variables!(rhs.x, allvars, maps[2])
        return MA.mutable_operate!(op, p, rhs)
    end
    get1(i) = (p.a[i], p.x.Z[i])
    get2(i) = (MA.scaling_convert(eltype(p.a), MA.copy_if_mutable(q.a[i])), copy(q.x.Z[i]))
    function set(i, t::_NoVarTerm)
        p.a[i] = t[1]
        p.x.Z[i] = t[2]
    end
    function push(t::_NoVarTerm)
        push!(p.a, t[1])
        push!(p.x.Z, t[2])
    end
    compare_monomials(t::_NoVarTerm, j::Int) = samevars_grlex(q.x.Z[j], t[2])
    compare_monomials(i::Int, j::Int) = compare_monomials(get1(i), j)
    combine(i::Int, j::Int) = p.a[i] = MA.operate!(op, p.a[i], q.a[j])
    combine(t::_NoVarTerm, j::Int) = (MA.operate!(op, t[1], q.a[j]), t[2])
    function resize(n)
        resize!(p.a, n)
        resize!(p.x.Z, n)
    end
    # We can modify the coefficient since it's the result of `combine`.
    keep(t::_NoVarTerm) = !MA.iszero!(t[1])
    keep(i::Int) = !MA.iszero!(p.a[i])
    MP.polynomial_merge!(
        nterms(p), nterms(q), get1, get2, set, push,
        compare_monomials, combine, keep, resize
    )
    return p
end

Base.:(+)(x::TermPoly{C}, y::TermPoly{C}) where C = plusorminus(x, y, +)
Base.:(-)(x::TermPoly{C}, y::TermPoly{C}) where C = plusorminus(x, y, -)
Base.:(+)(x::TermPoly{C}, y::Union{Monomial,PolyVar}) where C = x + Term{C}(y)
Base.:(+)(x::Union{Monomial,PolyVar}, y::TermPoly{C}) where C = Term{C}(x) + y

Base.:(-)(x::TermPoly{T}, y::DMonomialLike) where T = x - Term{T}(y)
Base.:(-)(x::DMonomialLike, y::TermPoly{T}) where T = Term{T}(x) - y

Base.:(-)(p::Polynomial) = Polynomial(-p.a, p.x)

include("mult.jl")
