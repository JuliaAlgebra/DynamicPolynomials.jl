export MonomialVector

# Invariant: Always sorted and no zero vector
struct MonomialVector{C} <: AbstractVector{Monomial{C}}
    vars::Vector{PolyVar{C}}
    Z::Vector{Vector{Int}}

    function MonomialVector{C}(
        vars::Vector{PolyVar{C}},
        Z::Vector{Vector{Int}},
    ) where {C}
        @assert !C || issorted(vars, rev = true)
        @assert all(z -> length(z) == length(vars), Z)
        @assert issorted(Z, lt = grlex)
        return new(vars, Z)
    end
end
function MonomialVector(
    vars::Vector{PolyVar{C}},
    Z::Vector{Vector{Int}},
) where {C}
    return MonomialVector{C}(vars, Z)
end
MonomialVector{C}() where {C} = MonomialVector{C}(PolyVar{C}[], Vector{Int}[])

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::MonomialVector)
    v = zeros(Bool, nvariables(m))
    for z in m.Z
        v = [v[i] || z[i] > 0 for i in eachindex(v)]
    end
    return MonomialVector(_vars(m)[v], Vector{Int}[z[v] for z in m.Z])
end

function Base.hash(m::MonomialVector, u::UInt)
    cm = canonical(m)
    if length(cm.Z) == 0
        hash([], u)
    elseif length(cm.Z) == 1
        hash(Monomial(_vars(cm), cm.Z[1]), u)
    else
        hash(_vars(cm), hash(cm.Z, hash(u)))
    end
end

# vars copied as it may be modifed by `_add_variables!`
# `mutable_copy` recursively copies the vector or vector of integers.
function MA.mutable_copy(m::MV) where {MV<:MonomialVector}
    return MV(copy(m.vars), MA.mutable_copy(m.Z))
end
Base.copy(m::MonomialVector) = MA.mutable_copy(m)
function Base.getindex(x::MonomialVector, I::AbstractVector{Bool})
    return typeof(x)(x.vars, x.Z[I])
end
function Base.getindex(x::MonomialVector, I::AbstractVector{<:Integer})
    return typeof(x)(x.vars, x.Z[sort(I)])
end
Base.getindex(x::MonomialVector, i::Integer) = Monomial(x.vars, x.Z[i])
Base.getindex(x::MonomialVector, i::CartesianIndex{1}) = x[i.I[1]]

function Base.deleteat!(x::MonomialVector, i)
    deleteat!(x.Z, i)
    return x
end
function Base.pop!(x::MonomialVector)
    pop!(x.Z)
    return x
end

Base.firstindex(x::MonomialVector) = firstindex(x.Z)
Base.lastindex(x::MonomialVector) = lastindex(x.Z)
Base.size(x::MonomialVector) = (length(x),)
Base.length(x::MonomialVector) = length(x.Z)
Base.isempty(x::MonomialVector) = length(x) == 0
Base.iterate(x::MonomialVector) = isempty(x) ? nothing : (x[1], 1)
function Base.iterate(x::MonomialVector, state::Int)
    return state < length(x) ? (x[state+1], state + 1) : nothing
end

function MultivariatePolynomials.extdegree(x::MonomialVector)
    return isempty(x) ? (0, 0) : extrema(sum.(x.Z))
end
function MultivariatePolynomials.mindegree(x::MonomialVector)
    return isempty(x) ? 0 : minimum(sum.(x.Z))
end
function MultivariatePolynomials.maxdegree(x::MonomialVector)
    return isempty(x) ? 0 : maximum(sum.(x.Z))
end

_vars(m::Union{Monomial,MonomialVector}) = m.vars

# Recognize arrays of monomials of this module
# [x, y] -> Vector{PolyVar}
# [x, x*y] -> Vector{Monomial}
# [1, x] -> Vector{Term{Int}}
const Dmonomial_vectorElemNonConstant{C} = Union{PolyVar{C},Monomial{C},Term{C}}
# [1] -> Vector{Int}
const Dmonomial_vectorElem{C} = Union{Int,Dmonomial_vectorElemNonConstant{C}}
const Dmonomial_vector{C} = AbstractVector{<:Dmonomial_vectorElem{C}}

function MP.empty_monomial_vector(vars::AbstractVector{PolyVar{C}}) where {C}
    return MonomialVector{C}(vars, Vector{Int}[])
end
function MP.empty_monomial_vector(t::Dmonomial_vectorElemNonConstant)
    return empty_monomial_vector(_vars(t))
end
function MP.empty_monomial_vector(
    ::Type{<:Dmonomial_vectorElemNonConstant{C}},
) where {C}
    return MonomialVector{C}()
end

function fillZfordeg!(Z, n, deg, ::Type{Val{true}}, filter::Function, ::Int)
    z = zeros(Int, n)
    z[end] = deg
    while true
        if filter(z)
            push!(Z, z)
            z = copy(z)
        end
        if z[1] == deg
            break
        end
        i = findfirst(i -> !iszero(z[i]), n:-1:2)
        j = (n:-1:2)[i]
        p = z[j]
        z[j] = 0
        z[end] = p - 1
        z[j-1] += 1
    end
end
function fillZrec!(Z, z, i, n, deg, filter::Function)
    if deg == 0
        if filter(z)
            push!(Z, copy(z))
        end
    else
        for i in i:i+n-1
            z[i] += 1
            fillZrec!(Z, z, i, n, deg - 1, filter)
            z[i] -= 1
        end
    end
end
function fillZfordeg!(
    Z,
    n,
    deg,
    ::Type{Val{false}},
    filter::Function,
    maxdeg::Int,
)
    z = zeros(Int, maxdeg * n - maxdeg + 1)
    start = length(Z) + 1
    fillZrec!(Z, z, 1, n, deg, filter)
    return reverse!(view(Z, start:length(Z)))
end
# List exponents in decreasing Graded Lexicographic Order
function getZfordegs(
    n,
    degs::AbstractVector{Int},
    ::Type{Val{C}},
    filter::Function,
) where {C}
    Z = Vector{Vector{Int}}()
    # For non-commutative, lower degree need to create a vector of exponent as large as for the highest degree
    maxdeg = isempty(degs) ? 0 : maximum(degs)
    for deg in sort(degs)
        fillZfordeg!(Z, n, deg, Val{C}, filter, maxdeg)
    end
    @assert issorted(Z, lt = grlex)
    return Z
end

function MonomialVector(
    vars::Vector{PolyVar{true}},
    degs::AbstractVector{Int},
    filter::Function = x -> true,
)
    vars = unique!(sort(vars, rev = true))
    return MonomialVector{true}(
        vars,
        getZfordegs(
            length(vars),
            degs,
            Val{true},
            z -> filter(Monomial(vars, z)),
        ),
    )
end

function getvarsforlength(vars::Vector{PolyVar{false}}, len::Int)
    n = length(vars)
    return map(i -> vars[((i-1)%n)+1], 1:len)
end
function MonomialVector(
    vars::Vector{PolyVar{false}},
    degs::AbstractVector{Int},
    filter::Function = x -> true,
)
    vars = unique!(sort(vars, rev = true))
    Z = getZfordegs(
        length(vars),
        degs,
        Val{false},
        z -> filter(Monomial(getvarsforlength(vars, length(z)), z)),
    )
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    return MonomialVector{false}(v, Z)
end
function MonomialVector(
    vars::Vector{<:PolyVar},
    degs::Int,
    filter::Function = x -> true,
)
    return MonomialVector(vars, [degs], filter)
end

function MP.monomials(vars::AbstractVector{<:PolyVar}, args...)
    return MonomialVector(vars, args...)
end
function MP.monomials(vars::Tuple{Vararg{PolyVar}}, args...)
    return monomials([vars...], args...)
end

#function MP.monomials(vars::TupOrVec{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
#    Z = getZfordegs(length(vars), degs, true, z -> filter(Monomial(vars, z)))
#    [Monomial{true}(vars, z) for z in Z]
#end
#function MP.monomials(vars::TupOrVec{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
#    Z = getZfordegs(length(vars), degs, false, z -> filter(Monomial(vars, z)))
#    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
#    [Monomial{false}(v, z) for z in Z]
#end
#MP.monomials(vars::TupOrVec{PV}, degs::Int, filter::Function = x->true) where {PV<:PolyVar} = monomials(vars, [degs], filter)

function buildZvarsvec(::Type{PV}, X::Dmonomial_vector) where {PV<:PolyVar}
    varsvec = Vector{PV}[
        (isa(x, Dmonomial_vectorElemNonConstant) ? _vars(x) : PolyVar[]) for
        x in X
    ]
    allvars, maps = mergevars(varsvec)
    nvars = length(allvars)
    Z = [zeros(Int, nvars) for i in 1:length(X)]
    offset = 0
    for (i, x) in enumerate(X)
        if isa(x, PolyVar)
            @assert length(maps[i]) == 1
            z = [1]
        elseif isa(x, Monomial)
            z = x.z
        elseif isa(x, Term)
            z = x.x.z
        else
            @assert isa(x, Int)
            z = Int[]
        end
        Z[i][maps[i]] = z
    end
    return allvars, Z
end

MP.sort_monomial_vector(X::MonomialVector) = (1:length(X), X)
function _sort_monomial_vector(X::Dmonomial_vector{C}) where {C}
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    σ = sortperm(Z, lt = grlex)
    return allvars, Z, σ
end
function _removedups!(Z::Vector{Vector{Int}}, σ::Vector{Int})
    dups = findall(i -> Z[σ[i]] == Z[σ[i-1]], 2:length(σ))
    return deleteat!(σ, dups)
end
function MP.sort_monomial_vector(X::Dmonomial_vector{C}) where {C}
    if isempty(X)
        Int[], MonomialVector{C}()
    else
        allvars, Z, σ = _sort_monomial_vector(X)
        _removedups!(Z, σ)
        σ, MonomialVector{C}(allvars, Z[σ])
    end
end

function MonomialVector{C}(X::Dmonomial_vector{C}) where {C}
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    sort!(Z, lt = grlex)
    dups = findall(i -> Z[i] == Z[i-1], 2:length(Z))
    deleteat!(Z, dups)
    return MonomialVector{C}(allvars, Z)
end
function MonomialVector(X)
    return monomial_vector_type(X)(X)
end

function MP.monomial_vector_type(
    X::Union{
        Dmonomial_vectorElemNonConstant{C},
        Type{<:Dmonomial_vectorElemNonConstant{C}},
        Dmonomial_vector{C},
        Type{<:Dmonomial_vector{C}},
    },
) where {C}
    return MonomialVector{C}
end
function MP.monomial_vector(X::Dmonomial_vector)
    return MonomialVector(X)
end
MP.monomial_vector(a, mv::MonomialVector) = (a, mv)

function MP.merge_monomial_vectors(ms::Vector{MonomialVector{C}}) where {C}
    m = length(ms)
    I = ones(Int, length(ms))
    L = length.(ms)
    X = Vector{Monomial{C}}()
    while any(I .<= L)
        min_monomial = nothing
        for i in 1:m
            if I[i] <= L[i]
                x = ms[i][I[i]]
                if min_monomial === nothing || min_monomial > x
                    min_monomial = x
                end
            end
        end
        @assert min_monomial !== nothing
        # to ensure that max is no more a union
        min_monomial === nothing && return X
        push!(X, min_monomial)
        for i in 1:m
            if I[i] <= L[i] && min_monomial == ms[i][I[i]]
                I[i] += 1
            end
        end
    end
    # There is no duplicate by construction
    return MonomialVector{C}(buildZvarsvec(PolyVar{C}, X)...)
end

function _add_variables!(
    monos::MonomialVector{C},
    allvars::Vector{PolyVar{C}},
    map,
) where {C}
    Future.copy!(monos.vars, allvars)
    if !isempty(monos.Z)
        tmp = similar(first(monos.Z))
        for z in monos.Z
            Future.copy!(tmp, z)
            resize!(z, length(allvars))
            fill!(z, 0)
            z[map] = tmp
        end
    end
end
