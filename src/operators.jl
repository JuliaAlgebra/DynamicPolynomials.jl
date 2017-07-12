import Base.(^), Base.(+)

# In Base/intfuncs.jl, x^p returns zero(x) when p == 0
# Since one(PolyVar) and one(Monomial) do not return
# a PolyVar and a Monomial, this results in type instability
# Defining the specific methods solve this problem and also make
# them a lot faster
^{C}(x::PolyVar{C}, i::Int) = Monomial{C}([x], [i])
^(x::Monomial{true}, i::Int) = Monomial{true}(x.vars, i*x.z)

myminivect{T}(x::T, y::T) = [x, y]
function myminivect{S,T}(x::S, y::T)
    U = promote_type(S, T)
    [U(x), U(y)]
end

function (+){C}(x::Term{C}, y::Term{C})
    if x.x == y.x
        Polynomial{C}([x.α+y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(y.α,x.α), [y.x,x.x])
    end
end

function (-){C}(x::Term{C}, y::Term{C})
    if x.x == y.x
        Polynomial{C}([x.α-y.α], [x.x])
    elseif x.x > y.x
        Polynomial{C}(myminivect(x.α,-y.α), [x.x,y.x])
    else
        Polynomial{C}(myminivect(-y.α,x.α), [y.x,x.x])
    end
end

(+){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) + Term(y)
(-){S<:Union{PolyVar,Monomial},T<:Union{PolyVar,Monomial}}(x::S, y::T) = Term(x) - Term(y)


function plusorminus{C, S, T}(p::TermPoly{C, S}, q::TermPoly{C, T}, isplus)
    varsvec = [_vars(p), _vars(q)]
    allvars, maps = myunion(varsvec)
    nvars = length(allvars)
    U = promote_type(S, T)
    a = Vector{U}()
    Z = Vector{Vector{Int}}()
    i = j = 1
    while i <= length(p) || j <= length(q)
        z = zeros(Int, nvars)
        if j > length(q) || (i <= length(p) && p[i].x > q[j].x)
            t = p[i]
            z[maps[1]] = t.x.z
            α = t.α
            i += 1
        elseif i > length(p) || q[j].x > p[i].x
            t = q[j]
            z[maps[2]] = t.x.z
            α = isplus ? t.α : -t.α
            j += 1
        else
            t = p[i]
            z[maps[1]] = t.x.z
            α = t.α
            s = q[j]
            α += isplus ? s.α : -s.α
            i += 1
            j += 1
        end
        push!(a, α)
        push!(Z, z)
    end

    Polynomial(a, MonomialVector{C}(allvars, Z))
end


(+){C}(x::TermPoly{C}, y::TermPoly{C}) = plusorminus(x, y, true)
(-){C}(x::TermPoly{C}, y::TermPoly{C}) = plusorminus(x, y, false)
(+){C, T, S<:Union{Monomial,PolyVar}}(x::TermPoly{C, T}, y::S) = x + Term{C, T}(y)
(+){C, T, S<:Union{Monomial,PolyVar}}(x::S, y::TermPoly{C, T}) = Term{C, T}(x) + y

(-){S<:Union{Monomial,PolyVar},T}(x::TermPoly{T}, y::S) = x - Term{T}(y)
(-){S<:Union{Monomial,PolyVar},T}(x::S, y::TermPoly{T}) = Term{T}(x) - y

(-)(p::Polynomial) = Polynomial(-p.a, p.x)

include("mult.jl")
