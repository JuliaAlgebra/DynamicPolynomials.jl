export Polynomial

# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
struct Polynomial{C, T} <: AbstractPolynomial{T}
    a::Vector{T}
    x::MonomialVector{C}

    function Polynomial{C, T}(a::Vector{T}, x::MonomialVector{C}) where {C, T}
        if length(a) != length(x) throw(ArgumentError("There should be as many coefficient than monomials"))
        end
        zeroidx = Int[]
        for (i,α) in enumerate(a)
            if iszero(α)
                push!(zeroidx, i)
            end
        end
        if !isempty(zeroidx)
            isnz = ones(Bool, length(a))
            isnz[zeroidx] = false
            nzidx = find(isnz)
            a = a[nzidx]
            x = x[nzidx]
        end
        new{C, T}(a, x)
    end
end

iscomm(::Type{Polynomial{C, T}}) where {C, T} = C

Base.copy(p::Polynomial{C, T}) where {C, T} = Polynomial{C, T}(copy(p.a), copy(p.x))
Base.zero(::Type{Polynomial{C, T}}) where {C, T} = Polynomial(T[], MonomialVector{C}())
Base.one(::Type{Polynomial{C, T}}) where {C, T} = Polynomial([one(T)], MonomialVector{C}(PolyVar{C}[], [Int[]]))
Base.zero(p::Polynomial{C, T}) where {C, T} = Polynomial(T[], emptymonovec(_vars(p)))
Base.one(p::Polynomial{C, T}) where {C, T} = Polynomial([one(T)], MonomialVector(_vars(p), [zeros(Int, nvariables(p))]))

Polynomial{C, T}(a::AbstractVector, x::MonomialVector) where {C, T} = Polynomial{C, T}(Vector{T}(a), x)
Polynomial{C, T}(a::AbstractVector, X::DMonoVec) where {C, T} = Polynomial{C, T}(monovec(a, X)...)
Polynomial{C}(a::Vector{T}, x) where {C, T} = Polynomial{C, T}(a, x)
Polynomial(af::Union{Function, Vector}, x::DMonoVec{C}) where {C} = Polynomial{C}(af, x)

Polynomial{C}(α) where {C} = Polynomial(Term{C}(α))
Polynomial(x::Union{PolyVar{C}, Monomial{C}}) where {C} = Polynomial(Term{C}(x))
Polynomial{C}(p::Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}) where {C} = Polynomial(p)

Polynomial(p::Polynomial) = p
Polynomial(t::Term{C, T}) where {C, T} = Polynomial{C, T}([t.α], [t.x])
Base.convert(::Type{Polynomial{C, T}}, α) where {C, T} = Polynomial(Term{C, T}(α))
Base.convert(::Type{Polynomial{C, T}}, m::DMonomialLike{C}) where {C, T} = Polynomial(Term{C, T}(m))
Base.convert(::Type{Polynomial{C, T}}, t::Term{C}) where {C, T} = Polynomial{C, T}([T(t.α)], [t.x])
Base.convert(::Type{Polynomial{C, T}}, p::Polynomial{C, T}) where {C, T} = p
Base.convert(::Type{Polynomial{C, T}}, p::Polynomial{C, S}) where {C, S, T} = Polynomial{C}(Vector{T}(p.a), p.x)
Base.convert(::Type{Polynomial{C, T}}, p::AbstractPolynomialLike) where {C, T} = Polynomial{C, T}(polynomial(p, T))

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

# FIXME why did I need it ?
Base.convert(::Type{Any}, p::Polynomial) = p

#Base.convert(::Type{PolyType{C}}, p::TermContainer{C}) where {C} = p

# needed to build [p Q; Q p] where p is a polynomial and Q is a matpolynomial in Julia v0.5
#Base.convert(::Type{TermType{C}}, p::TermContainer{C}) where {C} = p
#Base.convert(::Type{TermType{C, T}}, p::TermContainer{C, T}) where {C, T} = p

Base.endof(p::Polynomial) = length(p)
Base.length(p::Polynomial) = length(p.a)
Base.isempty(p::Polynomial) = isempty(p.a)
Base.start(::Polynomial) = 1
Base.done(p::Polynomial, state) = length(p) < state
Base.next(p::Polynomial, state) = (p[state], state+1)
#eltype(::Type{Polynomial{C, T}}) where {C, T} = T
Base.getindex(p::Polynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])

#Base.transpose(p::Polynomial) = Polynomial(map(transpose, p.a), p.x) # FIXME invalid age range update

struct TermIterator{C, T} <: AbstractVector{Term{C, T}}
    p::Polynomial{C, T}
end
Base.endof(p::TermIterator) = length(p.p)
Base.length(p::TermIterator) = length(p.p.a)
Base.size(p::TermIterator) = (length(p),)
Base.isempty(p::TermIterator) = isempty(p.p.a)
Base.start(::TermIterator) = 1
Base.done(p::TermIterator, state) = length(p.p) < state
Base.next(p::TermIterator, state) = (p.p[state], state+1)
Base.getindex(p::TermIterator, I::Int) = Term(p.p.a[I[1]], p.p.x[I[1]])

MP.terms(p::Polynomial) = TermIterator(p)
MP.coefficients(p::Polynomial) = p.a
MP.monomials(p::Polynomial) = p.x
_vars(p::Polynomial) = _vars(p.x)

MP.extdegree(p::Polynomial) = extdegree(p.x)
MP.mindegree(p::Polynomial) = mindegree(p.x)
MP.maxdegree(p::Polynomial) = maxdegree(p.x)

MP.leadingcoefficient(p::Polynomial) = first(p.a)
MP.leadingmonomial(p::Polynomial) = first(p.x)
MP.leadingterm(p::Polynomial) = first(p)

function MP.removeleadingterm(p::Polynomial)
    Polynomial(p.a[2:end], p.x[2:end])
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

function removedups(adup::Vector{T}, Zdup::Vector{Vector{Int}}) where {T}
    σ = sortperm(Zdup, rev=true, lt=grlex)
    Z = Vector{Vector{Int}}()
    a = Vector{T}()
    i = 0
    j = 1
    while j <= length(adup)
        k = σ[j]
        if j == 1 || Zdup[k] != Zdup[σ[j-1]]
            push!(Z, Zdup[k])
            push!(a, adup[k])
            i += 1
        else
            a[i] += adup[k]
        end
        j += 1
    end
    a, Z
end
function polynomialclean(vars::Vector{PolyVar{C}}, adup::Vector{T}, Zdup::Vector{Vector{Int}}) where {C, T}
    a, Z = removedups(adup, Zdup)
    Polynomial{C, T}(a, MonomialVector{C}(vars, Z))
end

MP.polynomial(a::AbstractVector, x::DMonoVec, s::MP.ListState) = Polynomial(a, x)

MP.polynomial(f::Function, x::AbstractVector) = Polynomial(f, x)
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
            Z = Vector{Vector{Int}}(N)
            a = Vector{T}(N)
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
            x = Vector{Monomial{C}}(N)
            a = Vector{T}(N)
            offset = 0
            for i in 1:n
                # for j in 1:n wouldn't be cache friendly for Q
                for j in i:n
                    k = trimap(i, j, n)
                    q = Q[i, j]
                    x[offset+k] = mv[i] * mv[j]
                    a[offset+k] = q
                    if i != j
                        offset += 1
                        x[offset+k] = mv[j] * mv[i]
                        a[offset+k] = q
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
