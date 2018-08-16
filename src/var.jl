export PolyVar, @polyvar, @ncpolyvar
export polyvecvar


function polyarrayvar(::Type{PV}, prefix, indices...) where {PV}
    map(i -> PV("$(prefix)" * join(i, "_")), Iterators.product(indices...))
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

struct PolyVar{C} <: AbstractVariable
    id::Int
    name::AbstractString
    function PolyVar{C}(name::AbstractString) where {C}
        # gensym returns something like Symbol("##42")
        # we first remove "##" and then parse it into an Int
        id = parse(Int, string(gensym())[3:end])
        new(id, name)
    end
end

Base.hash(x::PolyVar, u::UInt) = hash(x.id, u)
Base.broadcastable(x::PolyVar) = Ref(x)

MP.name(v::PolyVar) = v.name
MP.monomial(v::PolyVar) = Monomial(v)
_vars(v::PolyVar) = [v]

iscomm(::Type{PolyVar{C}}) where {C} = C

const VarVec{C} = Union{AbstractVector{PolyVar{C}}, NTuple{<:Integer, PolyVar{C}}}

function mergevars(varsvec::Vector{Vector{PV}}) where {PV<:PolyVar}
    n = length(varsvec)
    is = ones(Int, n)
    maps = zeros.(Int, length.(varsvec))
    nonempty = BitSet(findall(!isempty, varsvec))
    vars = Vector{PV}()
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
    vars, maps
end
