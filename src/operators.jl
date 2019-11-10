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

_getindex(p::Polynomial, i) = p[i]
_getindex(t::Term, i) = t
function plusorminus(p::TermPoly{C, S}, q::TermPoly{C, T}, op) where {C, S, T}
    varsvec = [_vars(p), _vars(q)]
    allvars, maps = mergevars(varsvec)
    nvars = length(allvars)
    U = Base.promote_op(op, S, T)
    a = Vector{U}()
    Z = Vector{Vector{Int}}()
    i = j = 1
    while i <= nterms(p) || j <= nterms(q)
        z = zeros(Int, nvars)
        if j > nterms(q) || (i <= nterms(p) && _getindex(p, i).x > _getindex(q, j).x)
            t = _getindex(p, i)
            z[maps[1]] = t.x.z
            α = convert(U, t.α)
            i += 1
        elseif i > nterms(p) || _getindex(q, j).x > _getindex(p, i).x
            t = _getindex(q, j)
            z[maps[2]] = t.x.z
            α = convert(U, op(t.α))
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

    Polynomial(a, MonomialVector{C}(allvars, Z))
end

function MA.mutable_operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial,
                             q::Union{PolyVar, Monomial, Term})
    return MA.mutable_operate!(op, p, polynomial(q))
end
function MA.mutable_operate!(op::Union{typeof(+), typeof(-)}, p::Polynomial,
                             q::Polynomial)
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
    i = j = 1
    # Invariant:
    # The terms p[0] -> p[i-1] are sorted and are smaller to the remaining terms.
    # The terms p[i] -> p[end] are sorted and still need to be added.
    # The terms q[j] -> p[end] are sorted and still need to be added.
    while i <= nterms(p) && j <= nterms(q)
        comp = samevars_grlex(q.x.Z[j], p.x.Z[i])
        if comp > 0
            break
        end
        if iszero(comp)
            p.a[i] = MA.operate!(op, p.a[i], q.a[j])
            j += 1
        end
        i += 1
    end
    if j <= nterms(q) && i <= nterms(p)
        # Invariant:
        # The terms in `tmp` are sorted and still need to be added.
        # Moreover, they are smaller than the terms p[i] -> p[end].
        tmp = Queue{Tuple{eltype(p.a), Vector{Int}}}()
        enqueue!(tmp, (p.a[i], p.x.Z[i]))
        p.a[i] = q.a[j]
        p.x.Z[i] = q.x.Z[j]
        i += 1
        j += 1
        while !isempty(tmp) && j <= nterms(q)
            α, z = front(tmp)
            comp = samevars_grlex(q.x.Z[j], z)
            if comp >= 0
                if comp > 0
                    α = q.a[j]
                    z = q.x.Z[j]
                else
                    α = MA.operate!(op, α, q.a[j])
                end
                j += 1
            end
            if comp <= 0
                dequeue!(tmp)
            end
            if i <= nterms(p)
                enqueue!(tmp, (p.a[i], p.x.Z[i]))
                p.a[i] = α
                p.x.Z[i] = z
            else
                push!(p.a, α)
                push!(p.x.Z, z)
            end
            i += 1
        end
        if !isempty(tmp)
            @assert j == nterms(q) + 1
            n = length(tmp)
            resize!(p.a, length(p.a) + n)
            resize!(p.x.Z, length(p.x.Z) + n)
            for k in i:(nterms(p) - n)
                p.a[k + n] = p.a[k]
                p.x.Z[k + n] = p.x.Z[k]
            end
            for k in i:(i + n - 1)
                p.a[k], p.x.Z[k] = dequeue!(tmp)
            end
        end
    end
    if j <= nterms(q)
        resize!(p.a, length(p.a) + nterms(q) - j + 1)
        resize!(p.x.Z, length(p.x.Z) + nterms(q) - j + 1)
        while j <= nterms(q)
            p.a[i] = q.a[j]
            p.x.Z[i] = q.x.Z[j]
            i += 1
            j += 1
        end
    end
    @assert j == nterms(q) + 1
    @assert i <= nterms(p) + 1
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
