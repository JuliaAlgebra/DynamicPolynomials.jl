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
    # We may assign a complex or real variable to its value (ordinary substitution).
    # We follow the following rules:
    # - If a single substitution rule determines the value of a real variable, just substitute it.
    # - If a single substitution rule determines the value of a complex variable or its conjugate, substitute the appropriate
    #   value whereever something related to this variable is found (i.e., the complex variable, the conjugate variable, or
    #   its real or imaginary part)
    # - If a single substitution rule determines the value of the real or imaginary part of a complex variable alone, then only
    #   replace the real or imaginary parts if they occur explicitly. Don't do a partial substitution, i.e., `z` with the rule
    #   `zᵣ => 1` is left alone and not changed into `1 + im*zᵢ`. Even if both the real and imaginary parts are substituted as
    #   two individual rules (which we don't know of in this method), `z` will not be replaced.
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
    # Every variable should be replaced by some value of type T
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
    # Some variable may not be replaced
    vals =
        Vector{promote_type(_substype(s), Variable{V,M})}(undef, length(vars))
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

_add_variables!(α, β) = α
_add_variables!(p::PolyType, α) = p
_add_variables!(α, p::PolyType) = MP.operate!!(*, α, one(p))
function _add_variables!(x::Variable, p::PolyType)
    return MP.operate!!(*, x, one(p))
end
function ___add_variables!(p, q)
    varsvec = [MP.variables(p), MP.variables(q)]
    allvars, maps = mergevars(varsvec)
    if length(allvars) != length(MP.variables(p))
        __add_variables!(p, allvars, maps[1])
    end
    return allvars, maps
end
function _add_variables!(p::PolyType, q::PolyType)
    ___add_variables!(p, q)
    return p
end

function _mono_eval(z::Vector{Int}, vals::AbstractVector)
    if length(z) != length(vals)
        error("Cannot evaluate a polynomial of `$(length(z))` variables with only `$(length(vals))` values.")
    end
    if isempty(z)
        return one(eltype(vals))^1
    end
    # `Base.power_by_squaring` does a `copy` if `z[1]` is `1`
    # which is redirected to `MA.mutable_copy`
    val = vals[1]^z[1]
    for i in 2:length(vals)
        if iszero(z[i])
            val = _add_variables!(val, vals[i])
        else
            val = MA.operate!!(*, val, vals[i]^z[i])
        end
    end
    return val
end

MP.substitute(::MP.AbstractSubstitutionType, ::Variable, vals::AbstractVector) = _mono_eval((1,), vals)
MP.substitute(::MP.AbstractSubstitutionType, m::Monomial, vals::AbstractVector) = _mono_eval(m.z, vals)
function MP.substitute(st::MP.AbstractSubstitutionType, t::_Term, vals::AbstractVector)
    return MP.coefficient(t) * MP.substitute(st, MP.monomial(t), vals)
end
function MP.substitute(
    ::MP.Eval,
    p::Polynomial{V,M,T},
    vals::AbstractVector{S},
) where {V,M,T,S}
    # I need to check for iszero otherwise I get : ArgumentError: reducing over an empty collection is not allowed
    if iszero(p)
        zero(MA.promote_operation(*, S, T))
    else
        sum(i -> p.a[i] * _mono_eval(p.x.Z[i], vals), eachindex(p.a))
    end
end
function MP.substitute(
    ::MP.Subs,
    p::Polynomial{V,M,T},
    vals::AbstractVector{S},
) where {V,M,T,S}
    Tout = MA.promote_operation(*, T, MP.coefficient_type(S))
    q = zero_with_variables(
        Polynomial{V,M,Tout},
        mergevars_of(Variable{V,M}, vals)[1],
    )
    for i in eachindex(p.a)
        MA.operate!(+, q, p.a[i] * _mono_eval(p.x.Z[i], vals))
    end
    return q
end

function MA.promote_operation(
    ::typeof(MP.substitute),
    ::Type{MP.Subs},
    ::Type{Monomial{V,M}},
    ::Type{Pair{Variable{V,M},T}},
) where {V,M,T}
    U = MA.promote_operation(*, T, T)
    return MA.promote_operation(*, U, Monomial{V,M})
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::PolyType,
    s::MP.AbstractSubstitution...,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), s))
end

# TODO resolve ambiguity. Can remove after:
# https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/pull/305
function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::PolyType,
    s::MP.AbstractMultiSubstitution,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), (s,)))
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::PolyType,
    s::MP.Substitutions,
)
    return MP.substitute(st, p, subsmap(st, MP.variables(p), s))
end

(v::Variable)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), v, s)
(m::Monomial)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), m, s)
(t::_Term)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), t, s)
(p::Polynomial)(s::MP.AbstractSubstitution...) = MP.substitute(MP.Eval(), p, s)

(p::Variable)(x::Number) = x
function (p::Monomial)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::Monomial)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, x)
end
(p::Monomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
function (p::_Term)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::_Term)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, x)
end
(p::_Term)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
function (p::Polynomial)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::Polynomial)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, x)
end
(p::Polynomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
