export Monomial

const TupOrVec{T} = Union{AbstractVector{T}, Tuple{Vararg{T}}}

# Invariant:
# vars is increasing
# z may contain 0's (otherwise, getindex of MonomialVector would be inefficient)
struct Monomial{C} <: AbstractMonomial
    vars::Vector{PolyVar{C}}
    z::Vector{Int}

    function Monomial{C}(vars::Vector{PolyVar{C}}, z::Vector{Int}) where {C}
        if length(vars) != length(z)
            throw(ArgumentError("There should be as many variables as exponents"))
        end
        new(vars, z)
    end
end

Monomial{C}(vars::Tuple{Vararg{PolyVar{C}}}, z::Vector{Int}) where C = Monomial{C}([vars...], z)

iscomm(::Type{Monomial{C}}) where C = C
Monomial{C}() where C = Monomial{C}(PolyVar{C}[], Int[])
Monomial(vars::TupOrVec{PolyVar{C}}, z::Vector{Int}) where C = Monomial{C}(vars, z)
function Base.convert(::Type{Monomial{C}}, x::PolyVar{C}) where C
    return Monomial{C}([x], [1])
end
Monomial(x::PolyVar{C}) where C = convert(Monomial{C}, x)
function MP.convertconstant(::Type{Monomial{C}}, α) where C
    α == 1 || error("Cannot convert $α to a Monomial{$C} as it is not one")
    Monomial{C}(PolyVar{C}[], Int[])
end

# defaults to commutative so that `Monomial(1)` is consistent with TypedPolynomials
Monomial(α::Number) = convert(Monomial{true}, α)

Base.broadcastable(m::Monomial) = Ref(m)
MA.mutable_copy(m::M) where {M<:Monomial} = M(copy(m.vars), copy(m.z))
Base.copy(m::Monomial) = MA.mutable_copy(m)

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::Monomial)
    list = m.z .> 0
    Monomial(_vars(m)[list], m.z[list])
end
function Base.hash(x::Monomial, u::UInt)
    cx = canonical(x)
    if iszero(nvariables(cx))
        hash(1, u)
    elseif nvariables(cx) == 1 && cx.z[1] == 1
        hash(cx.vars[1], u)
    else # TODO reduce power in MP
        hash(_vars(cx), hash(cx.z, u))
    end
end

MP.exponents(m::Monomial) = m.z
# /!\ vars not copied, do not mess with vars
_vars(m::Union{Monomial}) = m.vars

MP.monomial(m::Monomial) = m
# Does m1 divides m2 ?
#function MP.divides(m1::Monomial, m2::Monomial)
#    i = j = 1
#    while i <= length(m1.z) && j <= length(m2.z)
#        if m1.vars[i] == m2.vars[j]
#            if m1.z[i] > m2.z[j]
#                return false
#            end
#            i += 1
#            j += 1
#        elseif m1.vars[i] > m2.vars[j]
#            if !iszero(m1.z[i])
#                return false
#            end
#            i += 1
#        else
#            j += 1
#        end
#    end
#    i > length(m1.z)
#end
function _add_variables!(mono::Monomial, allvars, map)
    Future.copy!(mono.vars, allvars)
    tmp = copy(mono.z)
    resize!(mono.z, length(allvars))
    fill!(mono.z, 0)
    mono.z[map] = tmp
end
