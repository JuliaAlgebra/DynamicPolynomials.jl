# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
struct Polynomial{V,M,T} <: AbstractPolynomial{T}
    a::Vector{T}
    x::MonomialVector{V,M}

    function Polynomial{V,M,T}(
        a::Vector{T},
        x::MonomialVector{V,M},
    ) where {V,M,T}
        length(a) == length(x) || throw(
            ArgumentError("There should be as many coefficient than monomials"),
        )
        p = new{V,M,T}(a, x)
        _remove_zeros!(p)
        return p
    end
end
function Polynomial{V,M,T}(terms::AbstractVector{<:_Term{V,M}}) where {V,M,T}
    a = T[coefficient(t) for t in terms]
    monos = Monomial{V,M}[monomial(t) for t in terms]
    allvars, Z = buildZvarsvec(Variable{V,M}, monos)
    x = MonomialVector{V,M}(allvars, Z)
    return Polynomial{V,M,T}(a, x)
end

iscomm(::Type{Polynomial{V,M,T}}) where {V,M,T} = V, M

function _zero_with_variables(
    ::Type{Polynomial{V,M,T}},
    vars::Vector{Variable{V,M}},
) where {V,M,T}
    return Polynomial(T[], MonomialVector{V,M}(vars, Vector{Int}[]))
end

Base.broadcastable(p::Polynomial) = Ref(p)
function MA.mutable_copy(p::Polynomial{V,M,T}) where {V,M,T}
    return Polynomial{V,M,T}(MA.mutable_copy(p.a), MA.mutable_copy(p.x))
end
Base.copy(p::Polynomial) = MA.mutable_copy(p)
function Base.zero(::Type{Polynomial{V,M,T}}) where {V,M,T}
    return Polynomial(T[], MonomialVector{V,M}())
end
function Base.one(::Type{Polynomial{V,M,T}}) where {V,M,T}
    return Polynomial([one(T)], MonomialVector{V,M}(Variable{V,M}[], [Int[]]))
end
function Base.zero(p::Polynomial{V,M,T}) where {V,M,T}
    return Polynomial(T[], empty_monomial_vector(copy(MP.variables(p))))
end
function Base.one(p::Polynomial{V,M,T}) where {V,M,T}
    return Polynomial(
        [one(T)],
        MonomialVector(copy(MP.variables(p)), [zeros(Int, nvariables(p))]),
    )
end

function Polynomial{V,M,T}(a::AbstractVector, x::MonomialVector) where {V,M,T}
    return Polynomial{V,M,T}(Vector{T}(a), x)
end
function Polynomial{V,M,T}(a::AbstractVector, X::DMonoVec) where {V,M,T}
    return Polynomial{V,M,T}(monomial_vector(a, X)...)
end
Polynomial{V,M}(a::Vector{T}, x) where {V,M,T} = Polynomial{V,M,T}(a, x)
function Polynomial(af::Union{Function,Vector}, x::DMonoVec{V,M}) where {V,M}
    return Polynomial{V,M}(af, x)
end

# This is called by the default implementation of `Base.oneunit`, still in Julia v1.8 at least
Polynomial{V,M,T}(p::Polynomial{V,M,T}) where {V,M,T} = p

Base.convert(::Type{Polynomial{V,M,T}}, p::Polynomial{V,M,T}) where {V,M,T} = p
function Base.convert(
    ::Type{Polynomial{V,M,T}},
    t::AbstractTermLike,
) where {V,M,T}
    if iszero(t)
        _zero_with_variables(Polynomial{V,M,T}, variables(t))
    else
        # `exponents(::Variable)` gives a tuple
        return Polynomial{V,M,T}(
            [coefficient(t)],
            MonomialVector{V,M}(variables(t), [_vec(exponents(t))]),
        )
    end
end
function Base.convert(
    ::Type{Polynomial{V,M,T}},
    p::AbstractPolynomialLike,
) where {V,M,T}
    return Polynomial{V,M,T}(terms(p))
end

function Polynomial{V,M}(
    p::Union{Polynomial{V,M},_Term{V,M},Monomial{V,M},Variable{V,M}},
) where {V,M}
    return Polynomial(p)
end
Polynomial{V,M}(α) where {V,M} = Polynomial(_Term{V,M}(α))

Polynomial(p::Polynomial) = p
function Polynomial(t::_Term{V,M,T}) where {V,M,T}
    return convert(Polynomial{V,M,T}, mutable_copy(t))
end
function Polynomial(x::Union{Variable{V,M},Monomial{V,M}}) where {V,M}
    return Polynomial(_Term{V,M}(x))
end

#Base.convert(::Type{TermContainer{V,M,T}}, p::Polynomial{V,M}) where {V,M,T} = Polynomial{V,M,T}(p)

function Polynomial{V,M,T}(f::Function, x::MonomialVector{V,M}) where {V,M,T}
    a = T[f(i) for i in 1:length(x)]
    return Polynomial{V,M,T}(a, x)
end
function Polynomial{V,M,T}(f::Function, x::AbstractVector) where {V,M,T}
    σ, X = sort_monomial_vector(x)
    a = T[f(i) for i in σ]
    if length(x) > length(X)
        rev = Dict(X[j] => j for j in eachindex(σ))
        for i in eachindex(x)
            j = rev[x[i]]
            if i != σ[j]
                a[j] += f(i)
            end
        end
    end
    return Polynomial{V,M,T}(a, X)
end
function Polynomial{V,M}(f::Function, x) where {V,M}
    return Polynomial{V,M,Base.promote_op(f, Int)}(f, x)
end

#Base.convert(::Type{PolyType{V,M}}, p::TermContainer{V,M}) where {V,M} = p

# needed to build [p Q; Q p] where p is a polynomial and Q is a matpolynomial in Julia v0.5
#Base.convert(::Type{term_type{V,M}}, p::TermContainer{V,M}) where {V,M} = p
#Base.convert(::Type{term_type{V,M,T}}, p::TermContainer{V,M,T}) where {V,M,T} = p


#Base.transpose(p::Polynomial) = Polynomial(map(transpose, p.a), p.x) # FIXME invalid age range update

struct TermIterator{V,M,T} <: AbstractVector{_Term{V,M,T}}
    p::Polynomial{V,M,T}
end
Base.firstindex(p::TermIterator) = firstindex(p.p.a)
Base.lastindex(p::TermIterator) = lastindex(p.p.a)
Base.length(p::TermIterator) = length(p.p.a)
Base.size(p::TermIterator) = (length(p),)
Base.isempty(p::TermIterator) = isempty(p.p.a)
Base.iterate(p::TermIterator) = isempty(p) ? nothing : (p[1], 1)
function Base.iterate(p::TermIterator, state::Int)
    return state < length(p) ? (p[state+1], state + 1) : nothing
end

Base.getindex(p::TermIterator, I::Int) = MP.term(p.p.a[I[1]], p.p.x[I[1]])

MP.terms(p::Polynomial) = TermIterator(p)
MP.coefficients(p::Polynomial) = p.a
MP.monomials(p::Polynomial) = p.x
MP.variables(p::Polynomial) = MP.variables(p.x)

MP.extdegree(p::Polynomial) = extdegree(p.x)
MP.mindegree(p::Polynomial) = mindegree(p.x)
MP.maxdegree(p::Polynomial) = maxdegree(p.x)

function MP.leading_coefficient(p::Polynomial{V,M,T}) where {V,M,T}
    return iszero(p) ? zero(T) : last(p.a)
end
function MP.leading_monomial(p::Polynomial)
    return iszero(p) ? constant_monomial(p) : last(p.x)
end
MP.leading_term(p::Polynomial) = iszero(p) ? zero_term(p) : last(terms(p))

function MP.remove_leading_term(p::Polynomial)
    return Polynomial(p.a[1:end-1], p.x[1:end-1])
end
function MA.operate!(::typeof(MP.remove_leading_term), p::Polynomial)
    pop!(p.a)
    pop!(p.x)
    return p
end
function MP.remove_monomials(p::Polynomial, x::MonomialVector)
    # use the fact that monomials are sorted to do this O(n) instead of O(n^2)
    j = 1
    I = Int[]
    for (i, t) in enumerate(MP.terms(p))
        while j <= length(x) && x[j] < MP.monomial(t)
            j += 1
        end
        if j > length(x) || x[j] != MP.monomial(t)
            push!(I, i)
        end
    end
    return Polynomial(p.a[I], p.x[I])
end
function MP.remove_monomials(p::Polynomial, x::Vector)
    return remove_monomials(p, MonomialVector(x))
end

function _remove_zeros!(p::Polynomial)
    zeroidx = Int[]
    for (i, α) in enumerate(p.a)
        if iszero(α)
            push!(zeroidx, i)
        end
    end
    if !isempty(zeroidx)
        deleteat!(p.a, zeroidx)
        deleteat!(p.x.Z, zeroidx)
    end
end

function removedups_to!(
    a::Vector{T},
    Z::Vector{Vector{Int}},
    adup::Vector{T},
    Zdup::Vector{Vector{Int}},
    ::Type{M},
) where {T,M}
    _isless = let M = M
        (a, b) -> MP.compare(a, b, M) < 0
    end
    σ = sortperm(Zdup, lt = _isless)
    i = 0
    j = 1
    while j <= length(adup)
        k = σ[j]
        if j == 1 || Zdup[k] != Zdup[σ[j-1]]
            push!(Z, Zdup[k])
            push!(a, adup[k])
            i += 1
        else
            a[i] = MA.operate!!(+, a[i], adup[k])
        end
        j += 1
    end
end
function polynomialclean(
    vars::Vector{Variable{V,M}},
    adup::Vector{T},
    Zdup::Vector{Vector{Int}},
) where {V,M,T}
    Z = Vector{Int}[]
    a = T[]
    removedups_to!(a, Z, adup, Zdup, M)
    return Polynomial{V,M,T}(a, MonomialVector{V,M}(vars, Z))
end
function polynomialclean_to!(
    p::Polynomial{V,M,T},
    vars::Vector{Variable{V,M}},
    adup::Vector{T},
    Zdup::Vector{Vector{Int}},
) where {V,M,T}
    Future.copy!(p.x.vars, vars)
    empty!(p.a)
    empty!(p.x.Z)
    removedups_to!(p.a, p.x.Z, adup, Zdup, M)
    _remove_zeros!(p)
    return p
end

function MP.polynomial!(a::Vector, x::DMonoVec, ::MP.ListState)
    return Polynomial(a, x)
end
function MP.polynomial(a::AbstractVector, x::DMonoVec, s::MP.ListState)
    return MP.polynomial!(collect(a), MA.mutable_copy(x), s)
end

#MP.polynomial(f::Function, x::AbstractVector) = Polynomial(f, x)
#MP.polynomial(ts::AbstractVector{Term{V,M,T}}) where {V,M,T} = Polynomial(coefficient.(ts), monomial.(ts)) # FIXME invalid age range update

# i < j
function trimap(i, j, n)
    return div(n * (n + 1), 2) - div((n - i + 1) * (n - i + 2), 2) + j - i + 1
end
function MP.polynomial(Q::AbstractMatrix{T}, mv::MonomialVector) where {T}
    return MP.polynomial(Q, mv, Base.promote_op(+, T, T))
end
function MP.polynomial(
    Q::AbstractMatrix,
    mv::MonomialVector{V,M},
    ::Type{T},
) where {V,M,T}
    if isempty(Q)
        zero(Polynomial{V,M,T})
    else
        n = length(mv)
        if V <: Commutative
            N = trimap(n, n, n)
            Z = Vector{Vector{Int}}(undef, N)
            a = Vector{T}(undef, N)
            for i in 1:n
                for j in i:n
                    k = trimap(i, j, n)
                    Z[k] = mv.Z[i] + mv.Z[j]
                    if i == j
                        a[k] = Q[i, j]
                    else
                        a[k] = Q[i, j] + Q[j, i]
                    end
                end
            end
            v = MP.variables(mv)
        else
            N = n^2
            x = Vector{Monomial{V,M}}(undef, N)
            a = Vector{T}(undef, N)
            offset = 0
            for i in 1:n
                # for j in 1:n wouldn't be cache friendly for Q
                for j in i:n
                    k = trimap(i, j, n)
                    x[offset+k] = mv[i] * mv[j]
                    a[offset+k] = Q[i, j]
                    if i != j
                        offset += 1
                        x[offset+k] = mv[j] * mv[i]
                        a[offset+k] = Q[j, i]
                    end
                end
            end
            a, X = monomial_vector(a, x)
            v = MP.variables(X)
            Z = X.Z
        end
        polynomialclean(v, a, Z)
    end
end

function MA.operate!(::typeof(zero), p::Polynomial)
    empty!(p.a)
    empty!(p.x.Z)
    return p
end
function MA.operate!(::typeof(one), p::Polynomial{V,M,T}) where {V,M,T}
    if isempty(p.a)
        push!(p.a, one(T))
        push!(p.x.Z, zeros(Int, length(p.x.vars)))
    else
        resize!(p.a, 1)
        p.a[1] = MA.operate!!(one, p.a[1])
        resize!(p.x.Z, 1)
        MA.operate!(zero, p.x.Z[1])
    end
    return p
end

function MP.map_coefficients(f::Function, p::Polynomial; nonzero = false)
    return Polynomial(map(f, p.a), MA.mutable_copy(p.x))
end

function MP.map_coefficients!(f::Function, p::Polynomial; nonzero = false)
    map!(f, p.a, p.a)
    if !nonzero
        _remove_zeros!(p)
    end
    return p
end

function MP.map_coefficients_to!(
    output::Polynomial,
    f::Function,
    t::MP.AbstractTermLike;
    nonzero = false,
)
    return MP.map_coefficients_to!(output, f, polynomial(t); nonzero = nonzero)
end
function MP.map_coefficients_to!(
    output::Polynomial,
    f::Function,
    p::Polynomial;
    nonzero = false,
)
    resize!(output.a, length(p.a))
    map!(f, output.a, p.a)
    Future.copy!(output.x.vars, p.x.vars)
    # TODO reuse the part of `Z` that is already in `output`.
    resize!(output.x.Z, length(p.x.Z))
    for i in eachindex(p.x.Z)
        output.x.Z[i] = copy(p.x.Z[i])
    end
    if !nonzero
        _remove_zeros!(output)
    end
    return output
end

function MP.map_exponents(f::Function, p::Polynomial, m::DMonomialLike)
    return Polynomial(MA.mutable_copy(p.a), MP.map_exponents(f, p.x, m))
end

function MP.map_exponents!(f::Function, p::Polynomial, m::DMonomialLike)
    MP.map_exponents!(f, p.x, m)
    return p
end

function __add_variables!(p::Polynomial, allvars, map)
    return __add_variables!(p.x, allvars, map)
end
