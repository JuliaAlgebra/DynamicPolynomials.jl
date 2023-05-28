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
    # We may also assign a complex value to its conjugate, or just the real or imaginary parts
    # Any combination of z, conj(z), real(z), imag(z), imag(conj(z)) can occur in either the polynomial or the substitution,
    # and we must handle all of them correctly.
    # Note: This may or may not work... Issues can arise if the substitutions contain the real and imaginary (or only one of
    # those) of a variable separately whenever vals is not of the correct type:
    # - Unless subs() originally had a polynomial-valued rhs, vals will be scalars, monomials, or terms. So when we try to
    #   assign a polynomial to its value (which is necessary, as the one-step substitution of only the real or only the
    #   imaginary part is incomplete), this is an impossible conversion.
    # - The coefficients in vals might not be complex-valued; but to assign only a part of the variable, we necessarily need to
    #   introduce an explicit imaginary coefficient to the value.
    # Currently, we don't do anything to catch these errors.
    if s.first.kind == cpNone
        for j in eachindex(vars)
            if vars[j] == s.first
                vals[j] = s.second
                C == Commutative && break
            end
        end
    else
        for j in eachindex(vars)
            if vars[j].variable_order.order.id ==
               s.first.variable_order.order.id
                if s.first.kind == cpFull || s.first.kind == cpConj
                    value = s.first.kind == cpConj ? conj(s.second) : s.second
                    if vars[j].kind == cpFull
                        vals[j] = value
                    elseif vars[j].kind == cpConj
                        vals[j] = conj(value)
                    elseif vars[j].kind == cpReal
                        vals[j] = real(value)
                    else
                        vals[j] = imag(value)
                    end
                elseif s.first.kind == cpReal
                    isreal(s.second) || error(
                        "Cannot assign a complex value to the real part of an expression",
                    )
                    value = real(s.second) # just to make sure the type is correct
                    if vars[j].kind == cpFull
                        vals[j] = value + im * imag(vals[j])
                    elseif vars[j].kind == cpConj
                        vals[j] = value - im * imag(vals[j])
                    elseif vars[j].kind == cpReal
                        vals[j] = value
                    end
                    # else we know the real part but use the imaginary part; do nothing
                else
                    @assert(s.first.kind == cpImag)
                    isreal(s.second) || error(
                        "Cannot assign a complex value to the imaginary part of an expression",
                    )
                    value = real(s.second) # just to make sure the type is correct
                    if vars[j].kind == cpFull
                        vals[j] = real(vals[j]) + im * value
                    elseif vars[j].kind == cpConj
                        vals[j] = real(vals[j]) - im * value
                    elseif vars[j].kind == cpImag
                        vals[j] = value
                    end
                    # else we know the imaginary part but use the real part; do nothing
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
    return val
end

_subs(st, ::Variable, vals) = monoeval([1], vals::AbstractVector)
_subs(st, m::Monomial, vals) = monoeval(m.z, vals::AbstractVector)
function _subs(st, t::_Term, vals)
    return MP.coefficient(t) * monoeval(MP.monomial(t).z, vals::AbstractVector)
end
function _subs(
    ::MP.Eval,
    p::Polynomial{V,M,T},
    vals::AbstractVector{S},
) where {V,M,T,S}
    # I need to check for iszero otherwise I get : ArgumentError: reducing over an empty collection is not allowed
    if iszero(p)
        zero(Base.promote_op(*, S, T))
    else
        sum(i -> p.a[i] * monoeval(p.x.Z[i], vals), 1:length(p))
    end
end
function _subs(
    ::MP.Subs,
    p::Polynomial{V,M,T},
    vals::AbstractVector{S},
) where {V,M,T,S}
    Tout = MA.promote_operation(*, T, MP.coefficient_type(S))
    q = zero_with_variables(
        Polynomial{V,M,Tout},
        mergevars_of(Variable{V,M}, vals)[1],
    )
    for i in 1:length(p.a)
        MA.operate!(+, q, p.a[i] * monoeval(p.x.Z[i], vals))
    end
    return q
end

function MA.promote_operation(
    ::typeof(MP.substitute),
    ::Type{MP.Subs},
    ::Type{Monomial{V,M}},
    ::Type{Pair{Variable{V,M},T}},
) where {V,M,T}
    U = MA.promote_operation(^, T, Int)
    return _Term{V,M,U}
end

function MP.substitute(
    st::MP.AbstractSubstitutionType,
    p::PolyType,
    s::MP.Substitutions,
)
    return _subs(st, p, subsmap(st, MP.variables(p), s))
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
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
(p::Monomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
function (p::_Term)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::_Term)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
(p::_Term)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
function (p::Polynomial)(x::NTuple{N,<:Number}) where {N}
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
function (p::Polynomial)(x::AbstractVector{<:Number})
    return MP.substitute(MP.Eval(), p, variables(p) => x)
end
(p::Polynomial)(x::Number...) = MP.substitute(MP.Eval(), p, variables(p) => x)
