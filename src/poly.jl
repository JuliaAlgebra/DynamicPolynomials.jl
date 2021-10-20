export Polynomial

# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
struct Polynomial{C, T} <: AbstractPolynomial{T}
    a::Vector{T}
    x::MonomialVector{C}

    function Polynomial{C, T}(a::Vector{T}, x::MonomialVector{C}) where {C, T}
        length(a) == length(x) || throw(ArgumentError("There should be as many coefficient than monomials"))
        p = new{C, T}(a, x)
        _remove_zeros!(p)
        return p
    end
end

iscomm(::Type{Polynomial{C, T}}) where {C, T} = C

Base.broadcastable(p::Polynomial) = Ref(p)
MA.mutable_copy(p::Polynomial{C, T}) where {C, T} = Polynomial{C, T}(MA.mutable_copy(p.a), MA.mutable_copy(p.x))
Base.copy(p::Polynomial) = MA.mutable_copy(p)
Base.zero(::Type{Polynomial{C, T}}) where {C, T} = Polynomial(T[], MonomialVector{C}())
Base.one(::Type{Polynomial{C, T}}) where {C, T} = Polynomial([one(T)], MonomialVector{C}(PolyVar{C}[], [Int[]]))
Base.zero(p::Polynomial{C, T}) where {C, T} = Polynomial(T[], emptymonovec(copy(_vars(p))))
Base.one(p::Polynomial{C, T}) where {C, T} = Polynomial([one(T)], MonomialVector(copy(_vars(p)), [zeros(Int, nvariables(p))]))

Polynomial{C, T}(a::AbstractVector, x::MonomialVector) where {C, T} = Polynomial{C, T}(Vector{T}(a), x)
Polynomial{C, T}(a::AbstractVector, X::DMonoVec) where {C, T} = Polynomial{C, T}(monovec(a, X)...)
Polynomial{C}(a::Vector{T}, x) where {C, T} = Polynomial{C, T}(a, x)
Polynomial(af::Union{Function, Vector}, x::DMonoVec{C}) where {C} = Polynomial{C}(af, x)

# TODO Remove with MP v0.2.8
Polynomial{C, T}(p::Polynomial{C, T}) where {C, T} = p

Base.convert(::Type{Polynomial{C, T}}, p::Polynomial{C, T}) where {C, T} = p
function Base.convert(::Type{Polynomial{C, T}},
                      p::Polynomial{C, S}) where {C, S, T}
    return Polynomial{C}(convert(Vector{T}, p.a), p.x)
end
#function convert(::Type{Polynomial{C, T}},
#                 p::AbstractPolynomialLike) where {C, T}
#    return convert(Polynomial{C, T}, polynomial(p, T))
#end
function Base.convert(::Type{Polynomial{C, T}}, t::Term{C}) where {C, T}
    return Polynomial{C, T}(T[t.α], [t.x])
end
function Base.convert(::Type{Polynomial{C, T}}, m::DMonomialLike{C}) where {C, T}
    return Polynomial(convert(Term{C, T}, m))
end
function MP.convertconstant(::Type{Polynomial{C, T}}, α) where {C, T}
    return Polynomial(convert(Term{C, T}, α))
end

Polynomial{C}(p::Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}) where {C} = Polynomial(p)
Polynomial{C}(α) where {C} = Polynomial(Term{C}(α))

Polynomial(p::Polynomial) = p
Polynomial(t::Term{C, T}) where {C, T} = Polynomial{C, T}([t.α], [t.x])
Polynomial(x::Union{PolyVar{C}, Monomial{C}}) where {C} = Polynomial(Term{C}(x))

#Base.convert(::Type{TermContainer{C, T}}, p::Polynomial{C}) where {C, T} = Polynomial{C, T}(p)

function Polynomial{C, T}(f::Function, x::MonomialVector{C}) where {C, T}
    a = T[f(i) for i in 1:length(x)]
    Polynomial{C, T}(a, x)
end
function Polynomial{C, T}(f::Function, x::AbstractVector) where {C, T}
    σ, X = sortmonovec(x)
    a = T[f(i) for i in σ]
    Polynomial{C, T}(a, X)
end
Polynomial{C}(f::Function, x) where {C} = Polynomial{C, Base.promote_op(f, Int)}(f, x)

#Base.convert(::Type{PolyType{C}}, p::TermContainer{C}) where {C} = p

# needed to build [p Q; Q p] where p is a polynomial and Q is a matpolynomial in Julia v0.5
#Base.convert(::Type{TermType{C}}, p::TermContainer{C}) where {C} = p
#Base.convert(::Type{TermType{C, T}}, p::TermContainer{C, T}) where {C, T} = p

Base.length(p::Polynomial) = length(p.a)
Base.isempty(p::Polynomial) = isempty(p.a)
Base.iterate(p::Polynomial) = isempty(p) ? nothing : (p[1], 1)
function Base.iterate(p::Polynomial, state::Int)
    state < length(p) ? (p[state+1], state+1) : nothing
end
#eltype(::Type{Polynomial{C, T}}) where {C, T} = T
Base.getindex(p::Polynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])

#Base.transpose(p::Polynomial) = Polynomial(map(transpose, p.a), p.x) # FIXME invalid age range update

struct TermIterator{C, T} <: AbstractVector{Term{C, T}}
    p::Polynomial{C, T}
end
Base.firstindex(p::TermIterator) = firstindex(p.p.a)
Base.lastindex(p::TermIterator) = lastindex(p.p.a)
Base.length(p::TermIterator) = length(p.p.a)
Base.size(p::TermIterator) = (length(p),)
Base.isempty(p::TermIterator) = isempty(p.p.a)
Base.iterate(p::TermIterator) = isempty(p) ? nothing : (p[1], 1)
function Base.iterate(p::TermIterator, state::Int)
    state < length(p) ? (p[state+1], state+1) : nothing
end

Base.getindex(p::TermIterator, I::Int) = Term(p.p.a[I[1]], p.p.x[I[1]])

MP.terms(p::Polynomial) = TermIterator(p)
MP.coefficients(p::Polynomial) = p.a
MP.monomials(p::Polynomial) = p.x
_vars(p::Polynomial) = _vars(p.x)

MP.extdegree(p::Polynomial) = extdegree(p.x)
MP.mindegree(p::Polynomial) = mindegree(p.x)
MP.maxdegree(p::Polynomial) = maxdegree(p.x)

MP.leadingcoefficient(p::Polynomial{C, T}) where {C, T} = iszero(p) ? zero(T) : first(p.a)
MP.leadingmonomial(p::Polynomial) = iszero(p) ? constantmonomial(p) : first(p.x)
MP.leadingterm(p::Polynomial) = iszero(p) ? zeroterm(p) : first(terms(p))

function MP.removeleadingterm(p::Polynomial)
    Polynomial(p.a[2:end], p.x[2:end])
end
function MA.operate!(::typeof(MP.removeleadingterm), p::Polynomial)
    deleteat!(p.a, 1)
    deleteat!(p.x, 1)
    return p
end
function MP.removemonomials(p::Polynomial, x::MonomialVector)
    # use the fact that monomials are sorted to do this O(n) instead of O(n^2)
    j = 1
    I = Int[]
    for (i,t) in enumerate(p)
        while j <= length(x) && x[j] > t.x
            j += 1
        end
        if j > length(x) || x[j] != t.x
            push!(I, i)
        end
    end
    Polynomial(p.a[I], p.x[I])
end
MP.removemonomials(p::Polynomial, x::Vector) = removemonomials(p, MonomialVector(x))

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

function removedups_to!(a::Vector{T}, Z::Vector{Vector{Int}},
                        adup::Vector{T}, Zdup::Vector{Vector{Int}}) where T
    σ = sortperm(Zdup, rev=true, lt=grlex)
    i = 0
    j = 1
    while j <= length(adup)
        k = σ[j]
        if j == 1 || Zdup[k] != Zdup[σ[j-1]]
            push!(Z, Zdup[k])
            push!(a, adup[k])
            i += 1
        else
            a[i] = MA.operate!(+, a[i], adup[k])
        end
        j += 1
    end
end
function polynomialclean(vars::Vector{PolyVar{C}}, adup::Vector{T}, Zdup::Vector{Vector{Int}}) where {C, T}
    Z = Vector{Int}[]
    a = T[]
    removedups_to!(a, Z, adup, Zdup)
    Polynomial{C, T}(a, MonomialVector{C}(vars, Z))
end
function polynomialclean_to!(p::Polynomial{C, T}, vars::Vector{PolyVar{C}}, adup::Vector{T}, Zdup::Vector{Vector{Int}}) where {C, T}
    Future.copy!(p.x.vars, vars)
    empty!(p.a)
    empty!(p.x.Z)
    removedups_to!(p.a, p.x.Z, adup, Zdup)
    _remove_zeros!(p)
    return p
end

MP.polynomial!(a::Vector, x::DMonoVec, ::MP.ListState) = Polynomial(a, x)
MP.polynomial(a::AbstractVector, x::DMonoVec, s::MP.ListState) = MP.polynomial!(collect(a), MA.mutable_copy(x), s)

#MP.polynomial(f::Function, x::AbstractVector) = Polynomial(f, x)
#MP.polynomial(ts::AbstractVector{Term{C, T}}) where {C, T} = Polynomial(coefficient.(ts), monomial.(ts)) # FIXME invalid age range update

# i < j
function trimap(i, j, n)
    div(n*(n+1), 2) - div((n-i+1)*(n-i+2), 2) + j-i+1
end
MP.polynomial(Q::AbstractMatrix{T}, mv::MonomialVector) where T = MP.polynomial(Q, mv, Base.promote_op(+, T, T))
function MP.polynomial(Q::AbstractMatrix, mv::MonomialVector{C}, ::Type{T}) where {C, T}
    if isempty(Q)
        zero(Polynomial{C, T})
    else
        n = length(mv)
        if C
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
            v = _vars(mv)
        else
            N = n^2
            x = Vector{Monomial{C}}(undef, N)
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
            a, X = monovec(a, x)
            v = _vars(X)
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
function MA.operate!(::typeof(one), p::Polynomial{C, T}) where {C, T}
    if isempty(p.a)
        push!(p.a, one(T))
        push!(p.x.Z, zeros(Int, length(p.x.vars)))
    else
        resize!(p.a, 1)
        MA.operate!(one, p.a[1])
        resize!(p.x.Z, 1)
        MA.operate!(zero, p.x.Z[1])
    end
    return p
end
