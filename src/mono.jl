export Monomial

const TupOrVec{T} = Union{AbstractVector{T},Tuple{Vararg{T}}}

# Invariant:
# vars is increasing
# z may contain 0's (otherwise, getindex of MonomialVector would be inefficient)
struct Monomial{V,M} <: AbstractMonomial
    vars::Vector{Variable{V,M}}
    z::Vector{Int}

    function Monomial{V,M}(
        vars::Vector{Variable{V,M}},
        z::Vector{Int},
    ) where {V,M}
        if length(vars) != length(z)
            throw(
                ArgumentError("There should be as many variables as exponents"),
            )
        end
        return new(vars, z)
    end
end

function Monomial{V,M}(
    vars::Tuple{Vararg{Variable{V,M}}},
    z::Vector{Int},
) where {V,M}
    return Monomial{V,M}([vars...], z)
end

iscomm(::Type{<:Monomial{V}}) where {V} = iscomm(V)
Monomial{V,M}() where {V,M} = Monomial{V,M}(Variable{V,M}[], Int[])
function Monomial(vars::TupOrVec{Variable{V,M}}, z::Vector{Int}) where {V,M}
    return Monomial{V,M}(vars, z)
end
function Base.convert(::Type{Monomial{V,M}}, x::Variable{V,M}) where {V,M}
    return Monomial{V,M}([x], [1])
end
Monomial(x::Variable{V,M}) where {V,M} = convert(Monomial{V,M}, x)
function MP.convert_constant(::Type{Monomial{V,M}}, α) where {V,M}
    isone(α) ||
        error("Cannot convert `$α` to a `Monomial{$V,$M}` as it is not one")
    return Monomial{V,M}(Variable{V,M}[], Int[])
end

# defaults to commutative so that `Monomial(1)` is consistent with TypedPolynomials
function Monomial(α::Number)
    return convert(Monomial{Commutative{CreationOrder},Graded{LexOrder}}, α)
end

Base.broadcastable(m::Monomial) = Ref(m)
MA.mutable_copy(m::M) where {M<:Monomial} = M(copy(m.vars), copy(m.z))
Base.copy(m::Monomial) = MA.mutable_copy(m)

function MA.operate!(::typeof(constant_monomial), mono::Monomial)
    MA.operate!(zero, mono.z)
    return mono
end

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::Monomial)
    list = m.z .> 0
    return Monomial(MP.variables(m)[list], m.z[list])
end
function Base.hash(x::Monomial, u::UInt)
    cx = canonical(x)
    if iszero(nvariables(cx))
        hash(1, u)
    elseif nvariables(cx) == 1 && cx.z[1] == 1
        hash(cx.vars[1], u)
    else # TODO reduce power in MP
        hash(MP.variables(cx), hash(cx.z, u))
    end
end

MP.exponents(m::Monomial) = m.z
# /!\ vars not copied, do not mess with vars
MP.variables(m::Union{Monomial}) = m.vars

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

# for efficiency reasons
function Base.conj(x::Monomial{true})
    cv = conj.(x.vars)
    perm = sortperm(cv, rev = true)
    return Monomial{true}(cv[perm], x.z[perm])
end
