struct SafeValues{T}
    undefined::BitSet
    values::Vector{T}
end
function Base.setindex!(sv::SafeValues, v, i::Int)
    delete!(sv.undefined, i)
    sv.values[i] = v
end

function fillmap!(vals, vars::Vector{PolyVar{false}}, s::MP.Substitution)
    for j in eachindex(vars)
        if vars[j] == s.first
            vals[j] = s.second
        end
    end
end
function fillmap!(vals, vars::Vector{PolyVar{true}}, s::MP.Substitution)
    j = findfirst(isequal(s.first), vars)
    if j !== nothing
        vals[j] = s.second
    end
end
function fillmap!(vals, vars, s::MP.AbstractMultiSubstitution)
    for (v, x) in zip(s.first, s.second)
        fillmap!(vals, vars, v => x)
    end
end

function fillmap!(vals, vars, s1::MP.AbstractSubstitution, s2::MP.AbstractSubstitution...)
    fillmap!(vals, vars, s1)
    fillmap!(vals, vars, s2...)
end

_eltype(::T) where {T} = T
_eltype(t::Tuple) = Base.promote_typeof(t...)
_eltype(::Tuple{Vararg{T}}) where {T} = T
_eltype(::AbstractVector{T}) where {T} = T
_substype(s::MP.AbstractSubstitution) = _eltype(s.second)
_substype(s1::MP.AbstractSubstitution, s2::MP.AbstractSubstitution...) = promote_type(_substype(s1), _substype(s2...))
_substype(s::MP.Substitutions) = _substype(s...)

function _subsmap(::MP.Eval, vars, s::MP.Substitutions)
    # Every variable should be replaced by some value of type T
    vals = SafeValues(BitSet(1:length(vars)), Vector{_substype(s)}(undef, length(vars)))
    fillmap!(vals, vars, s...)
    if !isempty(vals.undefined)
        throw(ArgumentError("Variable `$(vars[first(vals.undefined)])` was not assigned a value. Use `subs` to substitute only a subset of the variables."))
    end
    return vals.values
end
function _subsmap(::MP.Subs, vars::Vector{PolyVar{C}}, s::MP.Substitutions) where {C}
    # Some variable may not be replaced
    vals = Vector{promote_type(_substype(s), PolyVar{C})}(undef, length(vars))
    Future.copy!(vals, vars)
    fillmap!(vals, vars, s...)
    return vals
end

subsmap(st, vars, s::MP.Substitutions) = _subsmap(st, vars, s)
_vec(a::AbstractVector) = a
_vec(a::Tuple) = [a...]
function subsmap(st, vars, s::Tuple{MP.VectorMultiSubstitution})
    if vars === s[1].first || vars == s[1].first # shortcut, === happens when the user do p(variables(p) => ...)
        _vec(s[1].second)
    else
        _subsmap(st, vars, s)
    end
end

function monoeval(z::Vector{Int}, vals::AbstractVector)
    @assert length(z) == length(vals)
    if isempty(z)
        return one(eltype(vals))^1
    end
    val = vals[1]^z[1]
    for i in 2:length(vals)
        if z[i] > 0
            val *= vals[i]^z[i]
        end
    end
    val
end

_subs(st, ::PolyVar, vals) = monoeval([1], vals::AbstractVector)
_subs(st, m::Monomial, vals) = monoeval(m.z, vals::AbstractVector)
_subs(st, t::Term, vals) = t.α * monoeval(t.x.z, vals::AbstractVector)
function _subs(::MP.Eval, p::Polynomial{C, T}, vals::AbstractVector{S}) where {C, T, S}
    # I need to check for iszero otherwise I get : ArgumentError: reducing over an empty collection is not allowed
    if iszero(p)
        zero(Base.promote_op(*, S, T))
    else
        sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p))
    end
end
function _subs(::MP.Subs, p::Polynomial{C, T}, vals::AbstractVector{S}) where {C, T, S}
    Tout = MA.promote_operation(*, T, MP.coefficienttype(S))
    q = zero_with_variables(Polynomial{C,Tout}, variables(p))
    for i in 1:length(p.a)
        MA.mutable_operate!(+, q, p.a[i] * monoeval(p.x.Z[i], vals))
    end
    return q
end

function MP.substitute(st::MP.AbstractSubstitutionType, p::PolyType, s::MP.Substitutions)
    _subs(st, p, subsmap(st, _vars(p), s))
end

(v::PolyVar)(s::MP.AbstractSubstitution...)    = MP.substitute(MP.Eval(), v, s)
(m::Monomial)(s::MP.AbstractSubstitution...)   = MP.substitute(MP.Eval(), m, s)
(t::Term)(s::MP.AbstractSubstitution...)       = MP.substitute(MP.Eval(), t, s)
(p::Polynomial)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), p, s)

(p::PolyVar)(x::Number) = x
(p::Monomial)(x::NTuple{N, <:Number}) where N = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Monomial)(x::AbstractVector{<:Number}) = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Monomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Term)(x::NTuple{N, <:Number}) where N = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Term)(x::AbstractVector{<:Number}) = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Term)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Polynomial)(x::NTuple{N, <:Number}) where N = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Polynomial)(x::AbstractVector{<:Number}) = MP.substitute(MP.Eval(), p, variables(p)=>x)
(p::Polynomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p)=>x)
