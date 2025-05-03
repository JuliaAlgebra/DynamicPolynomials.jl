export MonomialVector

# Invariant: Always sorted and no zero vector
struct MonomialVector{V,M} <: AbstractVector{Monomial{V,M}}
    vars::Vector{Variable{V,M}}
    Z::Vector{Vector{Int}}

    function MonomialVector{V,M}(
        vars::Vector{Variable{V,M}},
        Z::Vector{Vector{Int}},
    ) where {V,M}
        @assert !iscomm(V) || issorted(vars, rev = true)
        @assert all(z -> length(z) == length(vars), Z)

        return new{V,M}(vars, Z)
    end
end
function MonomialVector(
    vars::Vector{Variable{V,M}},
    Z::Vector{Vector{Int}},
) where {V,M}
    return MonomialVector{V,M}(vars, Z)
end
function MonomialVector{V,M}() where {V,M}
    return MonomialVector{V,M}(Variable{V,M}[], Vector{Int}[])
end

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::MonomialVector)
    v = zeros(Bool, nvariables(m))
    for z in m.Z
        v = [v[i] || z[i] > 0 for i in eachindex(v)]
    end
    return MonomialVector(MP.variables(m)[v], Vector{Int}[z[v] for z in m.Z])
end

function Base.hash(m::MonomialVector, u::UInt)
    cm = canonical(m)
    if length(cm.Z) == 0
        hash([], u)
    elseif length(cm.Z) == 1
        hash(Monomial(MP.variables(cm), cm.Z[1]), u)
    else
        hash(MP.variables(cm), hash(cm.Z, hash(u)))
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
# Complex-valued degrees for monomial vectors
for (fun, call, def, ret) in [
    (:extdegree_complex, :extrema, (0, 0), :((min(v1, v2), max(v1, v2)))),
    (:mindegree_complex, :minimum, 0, :(min(v1, v2))),
    (:maxdegree_complex, :maximum, 0, :(max(v1, v2))),
]
    eval(quote
        function MultivariatePolynomials.$fun(x::MonomialVector)
            isempty(x) && return $def
            vars = variables(x)
            @assert(!any(isrealpart, vars) && !any(isimagpart, vars))
            grouping = isconj.(vars)
            v1 = $call(sum(z[grouping]) for z in x.Z)
            grouping = map(!, grouping)
            v2 = $call(sum(z[grouping]) for z in x.Z)
            return $ret
        end
    end)
end
# faster complex-related functions
Base.isreal(x::MonomialVector) = all(isreal, x.vars)
Base.conj(x::MonomialVector) = MonomialVector(conj.(x.vars), x.Z)

MP.variables(m::Union{Monomial,MonomialVector}) = m.vars

# Recognize arrays of monomials of this module
# [x, y] -> Vector{Variable}
# [x, x*y] -> Vector{Monomial}
# [1, x] -> Vector{Term{Int}}
const DMonoVecElemNonConstant{V,M} =
    Union{Variable{V,M},Monomial{V,M},_Term{V,M}}
# [1] -> Vector{Int}
const DMonoVecElem{V,M} = Union{Int,DMonoVecElemNonConstant{V,M}}
const DMonoVec{V,M} = AbstractVector{<:DMonoVecElem{V,M}}

function MP.empty_monomial_vector(
    vars::AbstractVector{Variable{V,M}},
) where {V,M}
    return MonomialVector{V,M}(vars, Vector{Int}[])
end
function MP.empty_monomial_vector(t::DMonoVecElemNonConstant)
    return empty_monomial_vector(MP.variables(t))
end
function MP.empty_monomial_vector(
    ::Type{<:DMonoVecElemNonConstant{V,M}},
) where {V,M}
    return MonomialVector{V,M}()
end

# TODO replace by MP function
function _error_for_negative_degree(deg)
    if deg < 0
        throw(
            ArgumentError(
                "The degree should be a nonnegative number but the provided degree `$deg` is negative.",
            ),
        )
    end
end

function _fill_noncomm_exponents_rec!(Z, z, i, n, deg, ::Type{MP.LexOrder}, filter::Function)
    if deg == 0
        if filter(z)
            push!(Z, copy(z))
        end
    else
        for i in i:i+n-1
            z[i] += 1
            _fill_noncomm_exponents_rec!(Z, z, i, n, deg - 1, LexOrder, filter)
            z[i] -= 1
        end
    end
end

function _fill_exponents!(
    Z,
    n,
    deg,
    ::Type{NonCommutative},
    ::Type{MP.LexOrder},
    filter::Function,
    maxdeg::Int,
)
    _error_for_negative_degree(deg)
    _error_for_negative_degree(maxdeg)
    z = zeros(Int, maxdeg * n - maxdeg + 1)
    start = length(Z) + 1
    _fill_noncomm_exponents_rec!(Z, z, 1, n, deg, MP.LexOrder, filter)
    return reverse!(view(Z, start:length(Z)))
end

# List exponents in decreasing Graded Lexicographic Order
function getvarsforlength(vars::Vector{<:Variable{<:NonCommutative}}, len::Int)
    n = length(vars)
    return map(i -> vars[((i-1)%n)+1], 1:len)
end

function MonomialVector(
    vars::Vector{<:Variable{<:NonCommutative,M}},
    degs::AbstractVector{Int},
    filter::Function = x -> true,
) where {M}
    vars = unique!(sort(vars, rev = true))
    Z = _all_exponents(
        length(vars),
        degs,
        NonCommutative,
        M,
        z -> filter(Monomial(getvarsforlength(vars, length(z)), z)),
    )
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    return MonomialVector(v, Z)
end

function MonomialVector(
    vars::Vector{<:Variable{<:Commutative,M}},
    degs::AbstractVector{Int},
    filter::Function = x -> true,
) where {M}
    vars = unique!(sort(vars, rev = true))
    Z = Iterators.Filter(MP.ExponentsIterator{M}(
        zeros(Int, length(vars));
        mindegree = minimum(degs),
        maxdegree = maximum(degs),
    )) do z
        filter(Monomial(vars, z))
    end
    return MonomialVector(vars, collect(Z))
end

function MonomialVector(
    vars::Vector{<:Variable},
    degs::Int,
    filter::Function = x -> true,
)
    return MonomialVector(vars, [degs], filter)
end

function MP.monomials(vars::AbstractVector{<:Variable}, args...)
    return MonomialVector(vars, args...)
end

function MP.monomials(vars::Tuple{Vararg{Variable}}, args...)
    return monomials([vars...], args...)
end

#function MP.monomials(vars::TupOrVec{Variable{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
#    Z = _all_exponents(length(vars), degs, true, z -> filter(Monomial(vars, z)))
#    [Monomial{true}(vars, z) for z in Z]
#end
#function MP.monomials(vars::TupOrVec{<:Variable{<:NonCommutative}}, degs::AbstractVector{Int}, filter::Function = x->true)
#    Z = _all_exponents(length(vars), degs, false, z -> filter(Monomial(vars, z)))
#    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
#    [Monomial(v, z) for z in Z]
#end
#MP.monomials(vars::TupOrVec{PV}, degs::Int, filter::Function = x->true) where {PV<:Variable} = monomials(vars, [degs], filter)

function buildZvarsvec(::Type{PV}, X::DMonoVec) where {PV<:Variable}
    varsvec = Vector{PV}[
        (isa(x, DMonoVecElemNonConstant) ? MP.variables(x) : Variable[]) for
        x in X
    ]
    allvars, maps = mergevars(varsvec)
    nvars = length(allvars)
    Z = [zeros(Int, nvars) for i in 1:length(X)]
    for (i, x) in enumerate(X)
        if isa(x, Variable)
            @assert length(maps[i]) == 1
            z = [1]
        elseif isa(x, Monomial)
            z = x.z
        elseif isa(x, _Term)
            z = MP.monomial(x).z
        else
            @assert isa(x, Int)
            z = Int[]
        end
        Z[i][maps[i]] = z
    end
    return allvars, Z
end

MP.sort_monomial_vector(X::MonomialVector) = (1:length(X), X)
function _sort_monomial_vector(X::DMonoVec{V,M}) where {V,M}
    allvars, Z = buildZvarsvec(Variable{V,M}, X)
    _isless = let M = M
        (a, b) -> MP.compare(a, b, M) < 0
    end
    σ = sortperm(Z, lt = _isless)
    return allvars, Z, σ
end
function _removedups!(Z::Vector{Vector{Int}}, σ::Vector{Int})
    dups = findall(i -> Z[σ[i]] == Z[σ[i-1]], 2:length(σ))
    return deleteat!(σ, dups)
end
function MP.sort_monomial_vector(X::DMonoVec{V,M}) where {V,M}
    if isempty(X)
        Int[], MonomialVector{V,M}()
    else
        allvars, Z, σ = _sort_monomial_vector(X)
        _removedups!(Z, σ)
        σ, MonomialVector{V,M}(allvars, Z[σ])
    end
end

function MonomialVector{V,M}(X::DMonoVec{V,M}) where {V,M}
    allvars, Z = buildZvarsvec(Variable{V,M}, X)
    _isless = let M = M
        (a, b) -> cmp(M(), a, b) < 0
    end
    sort!(Z, lt = _isless)
    dups = findall(i -> Z[i] == Z[i-1], 2:length(Z))
    deleteat!(Z, dups)
    return MonomialVector{V,M}(allvars, Z)
end
function MonomialVector(X)
    return monomial_vector_type(X)(X)
end

function MP.monomial_vector_type(
    ::Union{
        DMonoVecElemNonConstant{V,M},
        Type{<:DMonoVecElemNonConstant{V,M}},
        DMonoVec{V,M},
        Type{<:DMonoVec{V,M}},
    },
) where {V,M}
    return MonomialVector{V,M}
end
function MP.monomial_vector(X::DMonoVec)
    return MonomialVector(X)
end
MP.monomial_vector(a, mv::MonomialVector) = (a, mv)

function MP.merge_monomial_vectors(ms::Vector{MonomialVector{V,M}}) where {V,M}
    m = length(ms)
    I = ones(Int, length(ms))
    L = length.(ms)
    X = Monomial{V,M}[]
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
    return MonomialVector{V,M}(buildZvarsvec(Variable{V,M}, X)...)
end

function __add_variables!(
    monos::MonomialVector{V,M},
    allvars::Vector{Variable{V,M}},
    map,
) where {V,M}
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
