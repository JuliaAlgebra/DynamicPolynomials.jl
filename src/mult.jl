# Monomial multiplication helper functions
function multiplyexistingvar(i::Int, z::Vector{Int})
    newz = copy(z)
    newz[i] += 1
    return newz
end
function multiplyexistingvar(i::Int, Z::Vector{Vector{Int}})
    return Vector{Int}[multiplyexistingvar(i, z) for z in Z]
end
function insertvar(
    v::Vector{Variable{V,M}},
    x::Variable{V,M},
    i::Int,
) where {V,M}
    n = length(v)
    I = 1:i-1
    J = i:n
    K = J .+ 1
    w = Vector{Variable{V,M}}(undef, n + 1)
    w[I] = v[I]
    w[i] = x
    w[K] = v[J]
    return w
end
function insertvar(z::Vector{Int}, x::Variable, i::Int)
    n = length(z)
    I = 1:i-1
    J = i:n
    K = J .+ 1
    newz = Vector{Int}(undef, n + 1)
    newz[I] = z[I]
    newz[i] = 1
    newz[K] = z[J]
    return newz
end
function insertvar(Z::Vector{Vector{Int}}, x::Variable, i::Int)
    return Vector{Int}[insertvar(z, x, i) for z in Z]
end

# Commutative and noncommutative monomial multiplication
include("cmult.jl")
include("ncmult.jl")

MP.left_constant_mult(α, x::Monomial) = MP.term(α, MA.mutable_copy(x))
