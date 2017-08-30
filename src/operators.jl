import Base.(^), Base.(+)

# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(PolyVar) and one(Monomial) do not return
# a PolyVar and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
^(x::PolyVar{C}, i::Int) where {C} = Monomial{C}([x], [i])
^(x::Monomial{true}, i::Int) = Monomial{true}(x.vars, i*x.z)

myminivect(x::T, y::T) where {T} = [x, y]
function myminivect(x::S, y::T) where {S,T}
    U = promote_type(S, T)
    [U(x), U(y)]
end

function (+)(x::Term{C}, y::Term{C}) where {C}
    if x.x == y.x
        Polynomial{C}([x.α+y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(y.α,x.α), [y.x,x.x])
    end
end

function (-)(x::Term{C}, y::Term{C}) where {C}
    if x.x == y.x
        Polynomial{C}([x.α-y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,-y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(-y.α,x.α), [y.x,x.x])
    end
end

(+)(x::DMonomialLike, y::DMonomialLike) = Term(x) + Term(y)
(-)(x::DMonomialLike, y::DMonomialLike) = Term(x) - Term(y)

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
            α = U(t.α)
            i += 1
        elseif i > nterms(p) || _getindex(q, j).x > _getindex(p, i).x
            t = _getindex(q, j)
            z[maps[2]] = t.x.z
            α = U(op(t.α))
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


(+)(x::TermPoly{C}, y::TermPoly{C}) where C = plusorminus(x, y, +)
(-)(x::TermPoly{C}, y::TermPoly{C}) where C = plusorminus(x, y, -)
(+)(x::TermPoly{C, T}, y::Union{Monomial,PolyVar}) where {C, T} = x + Term{C, T}(y)
(+)(x::Union{Monomial,PolyVar}, y::TermPoly{C, T}) where {C, T} = Term{C, T}(x) + y

(-)(x::TermPoly{T}, y::DMonomialLike) where T = x - Term{T}(y)
(-)(x::DMonomialLike, y::TermPoly{T}) where T = Term{T}(x) - y

(-)(p::Polynomial) = Polynomial(-p.a, p.x)

include("mult.jl")
