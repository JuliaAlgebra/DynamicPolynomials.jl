export MonomialVector

# Invariant: Always sorted and no zero vector
struct MonomialVector{C} <: AbstractVector{Monomial{C}}
    vars::Vector{PolyVar{C}}
    Z::Vector{Vector{Int}}

    function MonomialVector{C}(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}}) where {C}
        for z in Z
            if length(vars) != length(z)
                throw(ArgumentError("There should be as many vars than exponents"))
            end
        end
        @assert issorted(Z, rev=true, lt=grlex)
        new(vars, Z)
    end
end
MonomialVector(vars::Vector{PolyVar{C}}, Z::Vector{Vector{Int}}) where {C} = MonomialVector{C}(vars, Z)
(::Type{MonomialVector{C}})() where {C} = MonomialVector{C}(PolyVar{C}[], Vector{Int}[])
emptymonovec(vars::VarVec{C}) where {C} = MonomialVector{C}(vars, Vector{Int}[])

# Generate canonical reperesentation by removing variables that are not used
function canonical(m::MonomialVector)
    v = zeros(Bool, nvariables(m))
    for z in m.Z
        v = [v[i] || z[i] > 0 for i in eachindex(v)]
    end
    MonomialVector(_vars(m)[v], [z[v] for z in m.Z])
end

function Base.hash(m::MonomialVector, u::UInt)
    cm = canonical(m)
    if length(cm.Z) == 0
        hash(0, u)
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

Base.endof(x::MonomialVector) = length(x)
Base.size(x::MonomialVector) = (length(x),)
Base.length(x::MonomialVector) = length(x.Z)
Base.isempty(x::MonomialVector) = length(x) == 0
Base.start(::MonomialVector) = 1
Base.done(x::MonomialVector, state) = length(x) < state
Base.next(x::MonomialVector, state) = (x[state], state+1)

MultivariatePolynomials.extdeg(x::MonomialVector) = extrema(sum.(x.Z))
MultivariatePolynomials.mindeg(x::MonomialVector) = minimum(sum.(x.Z))
MultivariatePolynomials.maxdeg(x::MonomialVector) = maximum(sum.(x.Z))

_vars(m::Union{Monomial, MonomialVector}) = m.vars

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
function getZfordegs(n, degs::AbstractVector{Int}, C::Bool, filter::Function)
    Z = Vector{Vector{Int}}()
    for deg in sort(degs, rev=true)
        fillZfordeg!(Z, n, deg, Val{C}, filter)
    end
    @assert issorted(Z, rev=true, lt=grlex)
    Z
end
function MonomialVector(vars::Vector{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
    MonomialVector{true}(vars, getZfordegs(length(vars), degs, true, filter))
end
function getvarsforlength(vars::Vector{PolyVar{false}}, len::Int)
    n = length(vars)
    map(i -> vars[((i-1) % n) + 1], 1:len)
end

function MonomialVector(vars::Vector{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    MonomialVector{false}(v, Z)
end
MonomialVector(vars::Vector{PolyVar{C}}, degs::Int, filter::Function = x->true) where {C} = MonomialVector(vars, [degs], filter)

function MP.monomials(vars::TupOrVec{PolyVar{true}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, true, filter)
    [Monomial{true}(vars, z) for z in Z]
end
function MP.monomials(vars::TupOrVec{PolyVar{false}}, degs::AbstractVector{Int}, filter::Function = x->true)
    Z = getZfordegs(length(vars), degs, false, filter)
    v = isempty(Z) ? vars : getvarsforlength(vars, length(first(Z)))
    [Monomial{false}(v, z) for z in Z]
end
MP.monomials(vars::TupOrVec{PV}, degs::Int, filter::Function = x->true) where {PV<:PolyVar} = monomials(vars, [degs], filter)

# Recognize arrays of monomials of this module
# [x, y] -> Vector{PolyVar}
# [x, x*y] -> Vector{Monomial}
# [1, x] -> Vector{Term{Int}}
const DMonoVecElemNonConstant{C} = Union{PolyVar{C}, Monomial{C}, Term{C}}
# [1] -> Vector{Int}
const DMonoVecElem{C} = Union{Int, DMonoVecElemNonConstant{C}}
const DMonoVec{C} = AbstractVector{<:DMonoVecElem{C}}

function buildZvarsvec(::Type{PV}, X::DMonoVec) where {PV<:PolyVar}
    varsvec = Vector{PV}[ (isa(x, DMonoVecElemNonConstant) ? _vars(x) : PolyVar[]) for x in X ]
    allvars, maps = myunion(varsvec)
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
function MP.sortmonovec(X::DMonoVec{C}) where {C}
    if isempty(X)
        Int[], MonomialVector{C}()
    else
        allvars, Z = buildZvarsvec(PolyVar{C}, X)
        σ = sortperm(Z, rev=true, lt=grlex)
        σ, MonomialVector{C}(allvars, Z[σ])
    end
end

function MonomialVector{C}(X::DMonoVec{C}) where C
    allvars, Z = buildZvarsvec(PolyVar{C}, X)
    sort!(Z, rev=true, lt=grlex)
    MonomialVector{C}(allvars, Z)
end
function MonomialVector(X)
    monovectype(X)(X)
end

MP.monovectype(::Type{<:DMonoVecElemNonConstant{C}}) where {C} = MonomialVector{C}
MP.monovectype(X::DMonoVec{C}) where {C} = MonomialVector{C}
MP.monovec(::Type{<:DMonoVecElemNonConstant{C}}) where {C} = MonomialVector{C}()
function MP.monovec(X::DMonoVec)
    MonomialVector(X)
end
MP.monovec(a, mv::MonomialVector) = (a, mv)
function MP.monovec(a, x::DMonoVec)
    if length(a) != length(x)
        throw(ArgumentError("There should be as many coefficient than monomials"))
    end
    σ, X = sortmonovec(x)
    (a[σ], X)
end

function MP.mergemonovec(ms::Vector{MonomialVector{C}}) where {C}
    m = length(ms)
    I = ones(Int, length(ms))
    L = length.(ms)
    X = Vector{Monomial{C}}()
    while any(I .<= L)
        max = Nullable{Monomial{C}}()
        for i in 1:m
            if I[i] <= L[i]
                x = ms[i][I[i]]
                if isnull(max) || get(max) < x
                    max = Nullable(x)
                end
            end
        end
        @assert !isnull(max)
        push!(X, get(max))
        for i in 1:m
            if I[i] <= L[i] && get(max) == ms[i][I[i]]
                I[i] += 1
            end
        end
    end
    X
end
