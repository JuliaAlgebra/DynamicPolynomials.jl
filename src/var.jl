export PolyVar, @polyvar, @ncpolyvar, @polycvar
export polyvecvar


function polyarrayvar(::Type{PV}, prefix, indices...) where {PV}
    map(i -> PV("$(prefix)[$(join(i, ","))]"), Iterators.product(indices...))
end

function buildpolyvar(::Type{PV}, var) where {PV}
    if isa(var, Symbol)
        var, :($(esc(var)) = $PV($"$var"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
        (2 ≤ length(var.args)) || error("Expected $var to have at least one index set")
        varname = var.args[1]
        prefix = string(varname)
        varname, :($(esc(varname)) = polyarrayvar($PV, $prefix, $(esc.(var.args[2:end])...)))
    end
end

function buildpolyvars(::Type{PV}, args) where {PV}
    vars = Symbol[]
    exprs = []
    for arg in args
        var, expr = buildpolyvar(PV, arg)
        push!(vars, var)
        push!(exprs, expr)
    end
    vars, exprs
end

# Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
macro polyvar(args...)
    vars, exprs = buildpolyvars(PolyVar{true}, args)
    :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
end

macro ncpolyvar(args...)
    vars, exprs = buildpolyvars(PolyVar{false}, args)
    :($(foldl((x,y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
end

macro polycvar(args...)
    vars, exprs = buildpolyvars(PolyVarComplex{true}, args)
    return :($(foldl((x, y) -> :($x; $y), exprs, init=:())); $(Expr(:tuple, esc.(vars)...)))
end

@enum ComplexKind cpNone cpFull cpConj cpReal cpImag

struct PolyVar{C} <: AbstractVariable
    id::Int
    name::String
    kind::ComplexKind

    function PolyVar{C}(name::AbstractString, kind::ComplexKind=cpNone) where {C}
        # gensym returns something like Symbol("##42")
        # we first remove "##" and then parse it into an Int
        id = parse(Int, string(gensym())[3:end])
        new(id, convert(String, name), kind)
    end
    function PolyVar{C}(id::Int, name::String, kind::ComplexKind) where {C}
        new(id, name, kind)
    end
end

struct PolyVarComplex{C}
    PolyVarComplex{C}(name::AbstractString) where {C} = PolyVar{C}(name, cpFull)
end

Base.hash(x::PolyVar, u::UInt) = hash(x.id, u)
Base.broadcastable(x::PolyVar) = Ref(x)

MP.name(v::PolyVar) = v.name
function MP.name_base_indices(v::PolyVar)
    splits = split(v.name, r"[\[,\]]\s*", keepempty=false)
    if length(splits) == 1
        return v.name, Int[]
    else
        return splits[1], parse.(Int, splits[2:end])
    end
end

MP.monomial(v::PolyVar) = Monomial(v)
_vars(v::PolyVar) = [v]

iscomm(::Type{PolyVar{C}}) where {C} = C

MP.iscomplex(x::PolyVar{C}) where {C} = x.kind == cpFull || x.kind == cpConj
MP.isrealpart(x::PolyVar{C}) where {C} = x.kind == cpReal
MP.isimagpart(x::PolyVar{C}) where {C} = x.kind == cpImag
MP.isconj(x::PolyVar{C}) where {C} = x.kind == cpConj
MP.ordvar(x::PolyVar{C}) where {C} = x.kind == cpNone || x.kind == cpFull ? x : PolyVar{C}(x.id, x.name, cpFull)

function Base.conj(x::PolyVar{C}) where {C}
    if x.kind == cpFull
        return PolyVar{C}(x.id, x.name, cpConj)
    elseif x.kind == cpConj
        return PolyVar{C}(x.id, x.name, cpFull)
    else
        return x
    end
end
# for efficiency reasons
function Base.conj(x::Monomial{true})
    cv = conj.(x.vars)
    perm = sortperm(cv, rev=true)
    return Monomial{true}(cv[perm], x.z[perm])
end

function Base.real(x::PolyVar{C}) where {C}
    if x.kind == cpFull || x.kind == cpConj
        return PolyVar{C}(x.id, x.name, cpReal)
    else
        return x
    end
end
function Base.imag(x::PolyVar{C}) where {C}
    if x.kind == cpFull
        return PolyVar{C}(x.id, x.name, cpImag)
    elseif x.kind == cpConj
        return Term{C,Int}(-1, PolyVar{C}(x.id, x.name, cpImag))
    else
        return MA.Zero()
    end
end

function mergevars_to!(vars::Vector{PV}, varsvec::Vector{Vector{PV}}) where {PV<:PolyVar}
    empty!(vars)
    n = length(varsvec)
    is = ones(Int, n)
    maps = zeros.(Int, length.(varsvec))
    nonempty = BitSet(findall(!isempty, varsvec))
    while !isempty(nonempty)
        imin = 0
        for i in nonempty
            if imin == 0 || varsvec[i][is[i]] > varsvec[imin][is[imin]]
                imin = i
            end
        end
        var = varsvec[imin][is[imin]]
        push!(vars, var)
        for i in nonempty
            if var == varsvec[i][is[i]]
                maps[i][is[i]] = length(vars)
                if is[i] == length(varsvec[i])
                    pop!(nonempty, i)
                else
                    is[i] += 1
                end
            end
        end
    end
    return maps
end
function mergevars(varsvec::Vector{Vector{PV}}) where {PV<:PolyVar}
    vars = PV[]
    maps = mergevars_to!(vars, varsvec)
    return vars, maps
end
function mergevars_of(::Type{PolyVar{C}}, polys::AbstractVector) where {C}
    varsvec = Vector{PolyVar{C}}[variables(p) for p in polys if p isa PolyType]
    # TODO avoid computing `maps`
    return mergevars(varsvec)
end
