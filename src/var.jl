export PolyVar, @polyvar, @ncpolyvar
export polyvecvar

function polyvecvar(::Type{PV}, prefix, idxset) where {PV}
    [PV("$(prefix * string(i))") for i in idxset]
end

function buildpolyvar(::Type{PV}, var) where {PV}
    if isa(var, Symbol)
        :($(esc(var)) = $PV($"$var"))
    else
        isa(var, Expr) || error("Expected $var to be a variable name")
        Base.Meta.isexpr(var, :ref) || error("Expected $var to be of the form varname[idxset]")
        length(var.args) == 2 || error("Expected $var to have one index set")
        varname = var.args[1]
        prefix = string(var.args[1])
        idxset = esc(var.args[2])
        :($(esc(varname)) = polyvecvar($PV, $prefix, $idxset))
    end
end

# Variable vector x returned garanteed to be sorted so that if p is built with x then vars(p) == x
macro polyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(PolyVar{true}, arg) for arg in args])
end
macro ncpolyvar(args...)
    reduce((x,y) -> :($x; $y), :(), [buildpolyvar(PolyVar{false}, arg) for arg in args])
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
