export MonomialVector

# Invariant: Always sorted and no zero vector
struct MonomialVector{C} <: AbstractVector{Monomial{C}}
    vars::Vector{PolyVar{C}}
    Z::Vector{Vector{Int}}

    function MonomialVector{C}(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}}) where {C}
        @assert !C || issorted(vars, rev=true)
        @assert all(z -> length(z) == length(vars), Z)
        @assert issorted(Z, rev=true, lt=grlex)
        new(vars, Z)
    end
end
MonomialVector(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}}) where {C} = MonomialVector{C}(vars, Z)
MonomialVector{C}() where {C} = MonomialVector{C}(PolyVar{C}[], Vector{Int}[])

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::MonomialVector)
    v = zeros(Bool, nvariables(m))
    for z in m.Z
        v = [v[i] || z[i] > 0 for i in eachindex(v)]
    end
    MonomialVector(_vars(m)[v], Vector{Int}[z[v] for z in m.Z])
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

# /!\ vars not copied, do not mess with vars
Base.copy(m::MV) where {MV<:MonomialVector} = MV(m.vars, copy(m.Z))
function Base.getindex(x::MV, I) where {MV<:MonomialVector}
    MV(x.vars, x.Z[sort(I)])
end
Base.getindex(x::MonomialVector, i::Integer) = Monomial(x.vars, x.Z[i])

Base.firstindex(x::MonomialVector) = firstindex(x.Z)
Base.lastindex(x::MonomialVector) = lastindex(x.Z)
Base.size(x::MonomialVector) = (length(x),)
Base.length(x::MonomialVector) = length(x.Z)
Base.isempty(x::MonomialVector) = length(x) == 0
Base.iterate(x::MonomialVector) = isempty(x) ? nothing : (x[1], 1)
function Base.iterate(x::MonomialVector, state::Int)
    state < length(x) ? (x[state+1], state+1) : nothing
end

MultivariatePolynomials.extdegree(x::MonomialVector) = isempty(x) ? (0, 0) : extrema(sum.(x.Z))
MultivariatePolynomials.mindegree(x::MonomialVector) = isempty(x) ? 0 : minimum(sum.(x.Z))
MultivariatePolynomials.maxdegree(x::MonomialVector) = isempty(x) ? 0 : maximum(sum.(x.Z))

_vars(m::Union{Monomial, MonomialVector}) = m.vars

# Recognize arrays of monomials of this module
# [x, y] -> Vector{PolyVar}
# [x, x*y] -> Vector{Monomial}
# [1, x] -> Vector{Term{Int}}
const DMonoVecElemNonConstant{C} = Union{PolyVar{C}, Monomial{C}, Term{C}}
# [1] -> Vector{Int}
const DMonoVecElem{C} = Union{Int, DMonoVecElemNonConstant{C}}
const DMonoVec{C} = AbstractVector{<:DMonoVecElem{C}}

MP.emptymonovec(vars::VarVec{C}) where {C} = MonomialVector{C}(vars, Vector{Int}[])
MP.emptymonovec(t::DMonoVecElemNonConstant) = emptymonovec(_vars(t))
MP.emptymonovec(::Type{<:DMonoVecElemNonConstant{C}}) where {C} = MonomialVector{C}()

function fillZfordeg!(Z, n, deg, ::Type{Val{true}}, filter::Function)
    z = zeros(Int, n)
    z[1] = deg
    while true
        if filter(z)
            push!(Z, z)
            z = copy(z)
        end
        if z[end] == deg
            break
        end
        sum = 1
        for j in (n-1):-1:1
            if z[j] != 0
                z[j] -= 1
                z[j+1] += sum
                break
            else
                sum += z[j+1]
                z[j+1] = 0
            end
        end
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
            fillZrec!(Z, z, i, n, deg-1, filter)
            z[i] -= 1
        end
    end
end
function fillZfordeg!(Z, n, deg, ::Type{Val{false}}, filter::Function)
    z = zeros(Int, deg * n - deg + 1)
    fillZrec!(Z, z, 1, n, deg, filter)
end
# List exponents in decreasing Graded Lexicographic Order
function getZfordegs(n, degs::AbstractVector{Int}, ::Type{Val{C}}, filter::Function) where C
    Z = Vector{Vector{Int}}()
    for deg in sort(degs, rev=true)
        fillZfordeg!(Z, n, deg, Val{C}, filter)
    end
    @assert issorted(Z, rev=true, lt=grlex)
    Z
end

function MonomialVector(vars::Vector{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
    MonomialVector{true}(vars, getZfordegs(length(vars), degs, Val{true}, z -> filter(Monomial(vars, z))))
end

function getvarsforlength(vars::Vector{PolyVar{false}}, len::Int)
    n = length(vars)
    map(i -> vars[((i-1) % n) + 1], 1:len)
end
function MonomialVector(vars::Vector{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, Val{false}, z -> filter(Monomial(getvarsforlength(vars, length(z)), z)))
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    MonomialVector{false}(v, Z)
end
MonomialVector(vars::Vector{<:PolyVar}, degs::Int, filter::Function = x->true) = MonomialVector(vars, [degs], filter)

MP.monomials(vars::AbstractVector{<:PolyVar}, args...) = MonomialVector(vars, args...)
MP.monomials(vars::Tuple{Vararg{PolyVar}}, args...) = monomials([vars...], args...)

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

function buildZvarsvec(::Type{PV}, X::DMonoVec) where {PV<:PolyVar}
    varsvec = Vector{PV}[ (isa(x, DMonoVecElemNonConstant) ? _vars(x) : PolyVar[]) for x in X ]
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
    allvars, Z
end

MP.sortmonovec(X::MonomialVector) = (1:length(X), X)
function _sortmonovec(X::DMonoVec{C}) where {C}
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    σ = sortperm(Z, rev=true, lt=grlex)
    allvars, Z, σ
end
function _removedups!(Z::Vector{Vector{Int}}, σ::Vector{Int})
    dups = findall(i -> Z[σ[i]] == Z[σ[i-1]], 2:length(σ))
    deleteat!(σ, dups)
end
function MP.sortmonovec(X::DMonoVec{C}) where {C}
    if isempty(X)
        Int[], MonomialVector{C}()
    else
        allvars, Z, σ = _sortmonovec(X)
        _removedups!(Z, σ)
        σ, MonomialVector{C}(allvars, Z[σ])
    end
end

function MonomialVector{C}(X::DMonoVec{C}) where C
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    sort!(Z, rev=true, lt=grlex)
    dups = findall(i -> Z[i] == Z[i-1], 2:length(Z))
    deleteat!(Z, dups)
    MonomialVector{C}(allvars, Z)
end
function MonomialVector(X)
    monovectype(X)(X)
end

MP.monovectype(X::Union{DMonoVecElemNonConstant{C}, Type{<:DMonoVecElemNonConstant{C}}, DMonoVec{C}, Type{<:DMonoVec{C}}}) where {C} = MonomialVector{C}
function MP.monovec(X::DMonoVec)
    MonomialVector(X)
end
MP.monovec(a, mv::MonomialVector) = (a, mv)

function MP.mergemonovec(ms::Vector{MonomialVector{C}}) where {C}
    m = length(ms)
    I = ones(Int, length(ms))
    L = length.(ms)
    X = Vector{Monomial{C}}()
    while any(I .<= L)
        max = nothing
        for i in 1:m
            if I[i] <= L[i]
                x = ms[i][I[i]]
                if max === nothing || max < x
                    max = x
                end
            end
        end
        @assert max !== nothing
        # to ensure that max is no more a union
        max === nothing && return X
        push!(X, max)
        for i in 1:m
            if I[i] <= L[i] && max == ms[i][I[i]]
                I[i] += 1
            end
        end
    end
    # There is no duplicate by construction
    return MonomialVector{C}(buildZvarsvec(PolyVar{C}, X)...)
end
