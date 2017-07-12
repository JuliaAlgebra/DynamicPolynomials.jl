export Monomial

# Invariant:
# vars is increasing
# z may contain 0's (otherwise, getindex of MonomialVector would be inefficient)
type Monomial{C} <: AbstractMonomial
    vars::Vector{PolyVar{C}}
    z::Vector{Int}

    function Monomial{C}(vars::Vector{PolyVar{C}}, z::Vector{Int}) where {C}
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many vars than exponents"))
        end
        new(vars, z)
    end
end

iscomm{C}(::Type{Monomial{C}}) = C
(::Type{Monomial{C}}){C}() = Monomial{C}(PolyVar{C}[], Int[])
Monomial{C}(vars::Vector{PolyVar{C}}, z::Vector{Int}) = Monomial{C}(vars, z)
Monomial{C}(x::PolyVar{C}) = Monomial{C}(x)

Base.copy{M<:Monomial}(m::M) = M(m.vars, copy(m.z))
Base.convert{C}(::Type{Monomial{C}}, x::PolyVar{C}) = Monomial{C}([x], [1])

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::Monomial)
    list = m.z .> 0
    Monomial(_vars(m)[list], m.z[list])
end
function Base.hash(x::Monomial, u::UInt)
    cx = canonical(x)
    if nvars(cx) == 0
        hash(1, u)
    elseif nvars(cx) == 1 && cx.z[1] == 1
        hash(cx.vars[1], u)
    else
        hash(_vars(cx), hash(cx.z, u))
    end
end

# /!\ vars not copied, do not mess with vars
MP.exponents(m::Monomial) = m.z
_vars(m::Union{Monomial}) = m.vars

MP.monomial(m::Monomial) = m
# Does m1 divides m2 ?
function MP.divides(m1::Monomial, m2::Monomial)
    i = j = 1
    while i <= length(m1.z) && j <= length(m2.z)
        if m1.vars[i] == m2.vars[j]
            if m1.z[i] > m2.z[j]
                return false
            end
            i += 1
            j += 1
        elseif m1.vars[i] > m2.vars[j]
            if !iszero(m1.z[i])
                return false
            end
            i += 1
        else
            j += 1
        end
    end
    i > length(m1.z)
end
