function multiplyexistingvar(i::Int, z::Vector{Int})
    newz = copy(z)
    newz[i] += 1
    return newz
end
function multiplyexistingvar(i::Int, Z::Vector{Vector{Int}})
    return Vector{Int}[multiplyexistingvar(i, z) for z in Z]
end
function insertvar(
    v::Vector{Variable{V,M}},
    x::Variable{V,M},
    i::Int,
) where {V,M}
    n = length(v)
    I = 1:i-1
    J = i:n
    K = J .+ 1
    w = Vector{Variable{V,M}}(undef, n + 1)
    w[I] = v[I]
    w[i] = x
    w[K] = v[J]
    return w
end
function insertvar(z::Vector{Int}, x::Variable, i::Int)
    n = length(z)
    I = 1:i-1
    J = i:n
    K = J .+ 1
    newz = Vector{Int}(undef, n + 1)
    newz[I] = z[I]
    newz[i] = 1
    newz[K] = z[J]
    return newz
end
function insertvar(Z::Vector{Vector{Int}}, x::Variable, i::Int)
    return Vector{Int}[insertvar(z, x, i) for z in Z]
end

include("cmult.jl")
include("ncmult.jl")

MP.left_constant_mult(α, x::Monomial) = MP.term(α, MA.mutable_copy(x))

function zero_with_variables(
    ::Type{Polynomial{V,M,T}},
    vars::Vector{Variable{V,M}},
) where {V,M,T}
    return Polynomial(T[], empty_monomial_vector(vars))
end

# I do not want to cast x to TermContainer because that would force the promotion of eltype(q) with Int
function Base.:(*)(x::DMonomialLike, p::Polynomial)
    return Polynomial(MA.mutable_copy(p.a), x * p.x)
end
function Base.:(*)(x::DMonomialLike{<:NonCommutative}, p::Polynomial)
    # Order may change, e.g. y * (x + y) = y^2 + yx
    return Polynomial(
        monomial_vector(MA.mutable_copy(p.a), [x * m for m in p.x])...,
    )
end

function _term_poly_mult(
    t::_Term{V,M,S},
    p::Polynomial{V,M,T},
    op::Function,
) where {V,M,S,T}
    U = MA.promote_operation(op, S, T)
    if iszero(t)
        zero(Polynomial{V,M,U})
    else
        n = nterms(p)
        allvars, maps = mergevars([MP.monomial(t).vars, p.x.vars])
        nv = length(allvars)
        # Necessary to annotate the type in case it is empty
        Z = Vector{Int}[zeros(Int, nv) for i in 1:n]
        for i in 1:n
            Z[i][maps[1]] = MP.monomial(t).z
            Z[i][maps[2]] += p.x.Z[i]
        end
        Polynomial(op.(MP.coefficient(t), p.a), MonomialVector(allvars, Z))
    end
end
Base.:(*)(p::Polynomial, t::_Term) = _term_poly_mult(t, p, (α, β) -> β * α)
Base.:(*)(t::_Term, p::Polynomial) = _term_poly_mult(t, p, *)
_sumprod(a, b) = a * b + a * b
function _mul(
    ::Type{T},
    p::AbstractPolynomialLike,
    q::AbstractPolynomialLike,
) where {T}
    return _mul(T, polynomial(p), polynomial(q))
end
function _mul(
    ::Type{T},
    p::Polynomial{<:Commutative},
    q::Polynomial{<:Commutative},
) where {T}
    samevars = MP.variables(p) == MP.variables(q)
    if samevars
        allvars = copy(MP.variables(p))
    else
        allvars, maps = mergevars([MP.variables(p), MP.variables(q)])
    end
    N = length(p) * length(q)
    Z = Vector{Vector{Int}}(undef, N)
    a = Vector{T}(undef, N)
    i = 0
    for u in MP.terms(p)
        for v in MP.terms(q)
            if samevars
                z = MP.monomial(u).z + MP.monomial(v).z
            else
                z = zeros(Int, length(allvars))
                z[maps[1]] += MP.monomial(u).z
                z[maps[2]] += MP.monomial(v).z
            end
            i += 1
            Z[i] = z
            a[i] = MP.coefficient(u) * MP.coefficient(v)
        end
    end
    return allvars, a, Z
end
function Base.:(*)(
    p::Polynomial{V,M,S},
    q::Polynomial{V,M,T},
) where {V<:Commutative,M,S,T}
    PT = MA.promote_operation(*, typeof(p), typeof(q))
    if iszero(p) || iszero(q)
        zero(PT)
    else
        polynomialclean(_mul(MP.coefficient_type(PT), p, q)...)
    end
end
function MA.operate_to!(
    p::Polynomial{V,M,T},
    ::typeof(*),
    q1::MP.AbstractPolynomialLike,
    q2::MP.AbstractPolynomialLike,
) where {V<:NonCommutative,M,T}
    if iszero(q1) || iszero(q2)
        MA.operate!(zero, p)
    else
        ts = _Term{V,M,T}[]
        MP.mul_to_terms!(ts, q1, q2)
        # TODO do better than create tmp
        tmp = polynomial!(ts)
        Future.copy!(p.a, tmp.a)
        Future.copy!(p.x.vars, tmp.x.vars)
        Future.copy!(p.x.Z, tmp.x.Z)
        return p
    end
end
function MA.operate_to!(
    p::Polynomial{<:Commutative,M,T},
    ::typeof(*),
    q1::MP.AbstractPolynomialLike,
    q2::MP.AbstractPolynomialLike,
) where {M,T}
    if iszero(q1) || iszero(q2)
        MA.operate!(zero, p)
    else
        polynomialclean_to!(p, _mul(T, q1, q2)...)
    end
end
function MA.operate!(
    ::typeof(*),
    p::Polynomial{V,M},
    q::Polynomial{V,M},
) where {V,M}
    if iszero(q)
        return MA.operate!(zero, p)
    elseif nterms(q) == 1
        return MA.operate!(*, p, MP.leading_term(q))
    else
        return MA.operate_to!(p, *, p, q)
    end
end
