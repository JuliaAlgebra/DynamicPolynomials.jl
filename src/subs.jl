struct SafeValues{T}
    undefined::BitSet
    values::Vector{T}
end
function Base.setindex!(sv::SafeValues, v, i::Int)
    delete!(sv.undefined, i)
    return sv.values[i] = v
end

function fillmap!(
    vals,
    vars::Vector{<:Variable{C}},
    s::MP.Substitution,
) where {C}
    if s.first.kind == REAL
        for j in eachindex(vars)
            if vars[j] == s.first
                vals[j] = s.second
                C == Commutative && break
            end
        end
    else
        s.first.kind == COMPLEX || throw(
            ArgumentError(
                "Substitution with complex variables requires the ordinary_variable in the substitution specification",
            ),
        )
        for j in eachindex(vars)
            if vars[j].variable_order.order.id ==
               s.first.variable_order.order.id
                if vars[j].kind == COMPLEX
                    vals[j] = s.second
                elseif vars[j].kind == CONJ
                    vals[j] = conj(s.second)
                elseif vars[j].kind == REAL_PART
                    vals[j] = real(s.second)
                else
                    vals[j] = imag(s.second)
                end
            end
        end
    end
end
function fillmap!(vals, vars, s::MP.AbstractMultiSubstitution)
    for (v, x) in zip(s.first, s.second)
        fillmap!(vals, vars, v => x)
    end
end

function fillmap!(
    vals,
    vars,
    s1::MP.AbstractSubstitution,
    s2::MP.AbstractSubstitution...,
)
    fillmap!(vals, vars, s1)
    return fillmap!(vals, vars, s2...)
end

_eltype(::T) where {T} = T
_eltype(t::Tuple) = Base.promote_typeof(t...)
_eltype(::Tuple{Vararg{T}}) where {T} = T
_eltype(::AbstractVector{T}) where {T} = T
_substype(s::MP.AbstractSubstitution) = _eltype(s.second)
function _substype(s1::MP.AbstractSubstitution, s2::MP.AbstractSubstitution...)
    return promote_type(_substype(s1), _substype(s2...))
end
_substype(s::MP.Substitutions) = _substype(s...)

function _subsmap(::MP.Eval, vars, s::MP.Substitutions)
    vals = SafeValues(
        BitSet(1:length(vars)),
        Vector{_substype(s)}(undef, length(vars)),
    )
    fillmap!(vals, vars, s...)
    if !isempty(vals.undefined)
        throw(
            ArgumentError(
                "Variable `$(vars[first(vals.undefined)])` was not assigned a value. Use `subs` to substitute only a subset of the variables.",
            ),
        )
    end
    return vals.values
end
function _subsmap(
    ::MP.Subs,
    vars::Vector{Variable{V,M}},
    s::MP.Substitutions,
) where {V,M}
    vals =
        Vector{promote_type(_substype(s), Variable{V,M})}(undef, length(vars))
    copy!(vals, vars)
    fillmap!(vals, vars, s...)
    return vals
end

subsmap(st, vars, s::MP.Substitutions) = _subsmap(st, vars, s)
_vec(a::AbstractVector) = a
_vec(a::Tuple) = [a...]
function subsmap(st, vars, s::Tuple{MP.VectorMultiSubstitution})
    if vars === s[1].first || vars == s[1].first
        _vec(s[1].second)
    else
        _subsmap(st, vars, s)
    end
end

function _mono_eval(z::Union{Vector{Int},Tuple}, vals::AbstractVector)
    if length(z) != length(vals)
        error("Cannot evaluate a polynomial of `$(length(z))` variables with only `$(length(vals))` values.")
    end
    if isempty(z)
        return one(eltype(vals))^1
    end
    val = vals[1]^z[1]
    for i in 2:length(vals)
        if !iszero(z[i])
            val = MA.operate!!(*, val, vals[i]^z[i])
        end
    end
    return val
end

MP.substitute(::MP.AbstractSubstitutionType, ::Variable, vals::AbstractVector) = _mono_eval((1,), vals)
MP.substitute(::MP.AbstractSubstitutionType, m::Monomial, vals::AbstractVector) = _mono_eval(m.z, vals)

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::DMonomialLike,
    s::MP.AbstractSubstitution...,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), s))
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::DMonomialLike,
    s::MP.AbstractMultiSubstitution,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), (s,)))
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::DMonomialLike,
    s::MP.Substitutions,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), s))
end

(v::Variable)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), v, s)
(m::Monomial)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), m, s)

(p::Variable)(x::Number) = x
function (p::Monomial)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::Monomial)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, x)
end
(p::Monomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
