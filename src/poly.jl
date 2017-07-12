export Polynomial

# Invariant:
# a and x might be empty: meaning it is the zero polynomial
# a does not contain any zeros
# x is increasing in the monomial order (i.e. grlex)
type Polynomial{C, T} <: AbstractPolynomial{T}
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

iscomm{C, T}(::Type{Polynomial{C, T}}) = C

Base.copy{C, T}(p::Polynomial{C, T}) = Polynomial{C, T}(copy(p.a), copy(p.x))
Base.zero{C, T}(::Type{Polynomial{C, T}}) = Polynomial(T[], MonomialVector{C}())
Base.one{C, T}(::Type{Polynomial{C, T}}) = Polynomial([one(T)], MonomialVector{C}(PolyVar{C}[], [Int[]]))
Base.zero{C, T}(p::Polynomial{C, T}) = Polynomial(T[], emptymonovec(_vars(p)))
Base.one{C, T}(p::Polynomial{C, T}) = Polynomial([one(T)], MonomialVector(_vars(p), [zeros(Int, nvars(p))]))

(::Type{Polynomial{C, T}}){C, T}(a::AbstractVector, x::MonomialVector) = Polynomial{C, T}(Vector{T}(a), x)
(::Type{Polynomial{C, T}}){C, T}(a::AbstractVector, X::DMonoVec) = Polynomial{C, T}(monovec(a, X)...)
(::Type{Polynomial{C}}){C, T}(a::Vector{T}, x) = Polynomial{C, T}(a, x)
Polynomial{C}(af::Union{Function, Vector}, x::DMonoVec{C}) = Polynomial{C}(af, x)

(::Type{Polynomial{C}}){C}(α) = Polynomial(Term{C}(α))
Polynomial{C}(x::Union{PolyVar{C}, Monomial{C}}) = Polynomial(Term{C}(x))
(::Type{Polynomial{C}}){C}(p::Union{Polynomial{C}, Term{C}, Monomial{C}, PolyVar{C}}) = Polynomial(p)

Polynomial(p::Polynomial) = p
Polynomial{C, T}(t::Term{C, T}) = Polynomial{C, T}([t.α], [t.x])
Base.convert{C, T}(::Type{Polynomial{C, T}}, x) = Polynomial(Term{C, T}(x))
Base.convert{C, T}(::Type{Polynomial{C, T}}, t::Term{C}) = Polynomial{C, T}([T(t.α)], [t.x])
Base.convert{C, T}(::Type{Polynomial{C, T}}, p::Polynomial{C, T}) = p
Base.convert{C, S, T}(::Type{Polynomial{C, T}}, p::Polynomial{C, S}) = Polynomial{C}(Vector{T}(p.a), p.x)

#Base.convert{C, T}(::Type{TermContainer{C, T}}, p::Polynomial{C}) = Polynomial{C, T}(p)

function (::Type{Polynomial{C, T}}){C, T}(f::Function, x::MonomialVector{C})
    a = T[f(i) for i in 1:length(x)]
    Polynomial{C, T}(a, x)
end
function (::Type{Polynomial{C, T}}){C, T}(f::Function, x::Vector)
    σ, X = sortmonovec(x)
    a = T[f(i) for i in σ]
    Polynomial{C, T}(a, X)
end
(::Type{Polynomial{C}}){C}(f::Function, x) = Polynomial{C, Base.promote_op(f, Int)}(f, x)

# FIXME why did I need it ?
Base.convert(::Type{Any}, p::Polynomial) = p

#Base.convert{C}(::Type{PolyType{C}}, p::TermContainer{C}) = p

# needed to build [p Q; Q p] where p is a polynomial and Q is a matpolynomial in Julia v0.5
#Base.convert{C}(::Type{TermType{C}}, p::TermContainer{C}) = p
#Base.convert{C, T}(::Type{TermType{C, T}}, p::TermContainer{C, T}) = p

Base.endof(p::Polynomial) = length(p)
Base.length(p::Polynomial) = length(p.a)
Base.isempty(p::Polynomial) = isempty(p.a)
Base.start(::Polynomial) = 1
Base.done(p::Polynomial, state) = length(p) < state
Base.next(p::Polynomial, state) = (p[state], state+1)
#eltype{C, T}(::Type{Polynomial{C, T}}) = T
Base.getindex(p::Polynomial, I::Int) = Term(p.a[I[1]], p.x[I[1]])

MP.terms(p::Polynomial) = p
MP.monomials(p::Polynomial) = p.x
_vars(p::Polynomial) = _vars(p.x)

MP.extdeg(p::Polynomial) = extdeg(p.x)
MP.mindeg(p::Polynomial) = mindeg(p.x)
MP.maxdeg(p::Polynomial) = maxdeg(p.x)

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

function removedups{T}(adup::Vector{T}, Zdup::Vector{Vector{Int}})
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
function polynomialclean{C, T}(vars::Vector{PolyVar{C}}, adup::Vector{T}, Zdup::Vector{Vector{Int}})
    a, Z = removedups(adup, Zdup)
    Polynomial{C, T}(a, MonomialVector{C}(vars, Z))
end

MP.polynomial(a::AbstractVector, x::DMonoVec) = Polynomial(a, x)

function MP.polynomial{C, T}(Q::AbstractMatrix, mv::MonomialVector{C}, ::Type{T})
    if isempty(Q)
        zero(Polynomial{C, T})
    else
        n = length(mv)
        U = typeof(2*Q[1, 1] + Q[1, 1])
        if C
            N = MP.trimap(n, n, n)
            Z = Vector{Vector{Int}}(N)
            a = Vector{U}(N)
            for i in 1:n
                for j in i:n
                    k = MP.trimap(i, j, n)
                    Z[k] = mv.Z[i] + mv.Z[j]
                    if i == j
                        a[k] = Q[i, j]
                    else
                        a[k] = 2*Q[i, j]
                    end
                end
            end
            v = _vars(mv)
        else
            N = n^2
            x = Vector{Monomial{C}}(N)
            a = Vector{U}(N)
            offset = 0
            for i in 1:n
                # for j in 1:n wouldn't be cache friendly for Q
                for j in i:n
                    k = MP.trimap(i, j, n)
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
