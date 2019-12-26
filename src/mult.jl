function multiplyexistingvar(v::Vector{PolyVar{C}}, x::PolyVar{C}, i::Int) where {C}
    updatez = z -> begin
        newz = copy(z)
        newz[i] += 1
        newz
    end
    copy(v), updatez
end
function insertvar(v::Vector{PolyVar{C}}, x::PolyVar{C}, i::Int) where {C}
    n = length(v)
    I = 1:i-1
    J = i:n
    K = J.+1
    w = Vector{PolyVar{C}}(undef, n+1)
    w[I] = v[I]
    w[i] = x
    w[K] = v[J]
    updatez = z -> begin
        newz = Vector{Int}(undef, n+1)
        newz[I] = z[I]
        newz[i] = 1
        newz[K] = z[J]
        newz
    end
    w, updatez
end

include("cmult.jl")
include("ncmult.jl")

MP.multconstant(α, x::Monomial)   = Term(α, MA.mutable_copy(x))
MP.mapcoefficientsnz(f::Function, p::Polynomial) = Polynomial(map(f, p.a), MA.mutable_copy(p.x))
function MP.mapcoefficientsnz_to!(output::Polynomial, f::Function, t::MP.AbstractTermLike)
    MP.mapcoefficientsnz_to!(output, f, polynomial(t))
end
function MP.mapcoefficientsnz_to!(output::Polynomial, f::Function, p::Polynomial)
    resize!(output.a, length(p.a))
    @. output.a = f(p.a)
    Future.copy!(output.x.vars, p.x.vars)
    # TODO reuse the part of `Z` that is already in `output`.
    resize!(output.x.Z, length(p.x.Z))
    for i in eachindex(p.x.Z)
        output.x.Z[i] = copy(p.x.Z[i])
    end
    return output
end

# I do not want to cast x to TermContainer because that would force the promotion of eltype(q) with Int
function Base.:(*)(x::DMonomialLike, p::Polynomial)
    Polynomial(MA.mutable_copy(p.a), x*p.x)
end
function Base.:(*)(x::DMonomialLike{false}, p::Polynomial)
    # Order may change, e.g. y * (x + y) = y^2 + yx
    Polynomial(monovec(MA.mutable_copy(p.a), [x*m for m in p.x])...)
end
function Base.:(*)(p::Polynomial, x::DMonomialLike)
    Polynomial(MA.mutable_copy(p.a), p.x*x)
end

function _term_poly_mult(t::Term{C, S}, p::Polynomial{C, T}, op::Function) where {C, S, T}
    U = MA.promote_operation(op, S, T)
    if iszero(t)
        zero(Polynomial{C, U})
    else
        n = nterms(p)
        allvars, maps = mergevars([t.x.vars, p.x.vars])
        nv = length(allvars)
        # Necessary to annotate the type in case it is empty
        Z = Vector{Int}[zeros(Int, nv) for i in 1:n]
        for i in 1:n
            Z[i][maps[1]] = t.x.z
            Z[i][maps[2]] += p.x.Z[i]
        end
        Polynomial(op.(t.α, p.a), MonomialVector(allvars, Z))
    end
end
Base.:(*)(p::Polynomial, t::Term) = _term_poly_mult(t, p, (α, β) -> β * α)
Base.:(*)(t::Term, p::Polynomial) = _term_poly_mult(t, p, *)
_sumprod(a, b) = a * b + a * b
function _mul(::Type{T}, p::AbstractPolynomialLike, q::AbstractPolynomialLike) where T
    return _mul(T, polynomial(p), polynomial(q))
end
function _mul(::Type{T}, p::Polynomial{true}, q::Polynomial{true}) where T
    samevars = _vars(p) == _vars(q)
    if samevars
        allvars = copy(_vars(p))
    else
        allvars, maps = mergevars([_vars(p), _vars(q)])
    end
    N = length(p) * length(q)
    Z = Vector{Vector{Int}}(undef, N)
    a = Vector{T}(undef, N)
    i = 0
    for u in p
        for v in q
            if samevars
                z = u.x.z + v.x.z
            else
                z = zeros(Int, length(allvars))
                z[maps[1]] += u.x.z
                z[maps[2]] += v.x.z
            end
            i += 1
            Z[i] = z
            a[i] = u.α * v.α
        end
    end
    return allvars, a, Z
end
function Base.:(*)(p::Polynomial{true, S}, q::Polynomial{true, T}) where {S, T}
    if iszero(p) || iszero(q)
        zero(MA.promote_operation(*, typeof(p), typeof(q)))
    else
        polynomialclean(_mul(MA.promote_operation(MA.add_mul, S, T), p, q)...)
    end
end
function MA.mutable_operate_to!(p::Polynomial{false, T}, ::typeof(*), q1::MP.AbstractPolynomialLike, q2::MP.AbstractPolynomialLike) where T
    if iszero(q1) || iszero(q2)
        MA.mutable_operate!(zero, p)
    else
        ts = Term{false, T}[]
        MP.mul_to_terms!(ts, q1, q2)
        # TODO do better than create tmp
        tmp = polynomial!(ts)
        Future.copy!(p.a, tmp.a)
        Future.copy!(p.x.vars, tmp.x.vars)
        Future.copy!(p.x.Z, tmp.x.Z)
        return p
    end
end
function MA.mutable_operate_to!(p::Polynomial{true, T}, ::typeof(*), q1::MP.AbstractPolynomialLike, q2::MP.AbstractPolynomialLike) where T
    if iszero(q1) || iszero(q2)
        MA.mutable_operate!(zero, p)
    else
        polynomialclean_to!(p, _mul(T, q1, q2)...)
    end
end
function MA.mutable_operate!(::typeof(*), p::Polynomial{C}, q::Polynomial{C}) where C
    return MA.mutable_operate_to!(p, *, p, q)
end
