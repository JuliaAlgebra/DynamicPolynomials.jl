# Dynamic Polynomials

| **Build Status** | **References to cite** |
|:----------------:|:----------------------:|
| [![Build Status][build-img]][build-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![Codecov branch][codecov-img]][codecov-url] | |

Sparse dynamic representation of multivariate polynomials that can be used with [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) (see the documentation there for more information).
Both commutative and non-commutative variables are supported.
The following types are defined:

* `Variable{V,M}`: A variable which is commutative with `*` when `V<:Commutative`. Commutative variables are created using the `@polyvar` macro, e.g. `@polyvar x y`, `@polyvar x[1:8]` and non-commutative variables are created likewise using the `@ncpolyvar` macro. The type parameter `M` is the monomial ordering.
* `Monomial{V,M}`: A product of variables: e.g. `x*y^2`.
* `MultivariatePolynomials.Term{T,Monomial{V,M}}`: A product between an element of type `T` and a `Monomial{V,M}`, e.g `2x`, `3.0x*y^2`.
* `Polynomial{V,M,T}`: A sum of `Term{T,Monomial{V,M}}`, e.g. `2x + 3.0x*y^2 + y`.

All common algebraic operations between those types are designed to be as efficient as possible without doing any assumption on `T`.
Typically, one imagine `T` to be a subtype of `Number` but it can be anything.
This is useful for example in the package [PolyJuMP](https://github.com/jump-dev/PolyJuMP.jl) where `T` is often an affine expression of [JuMP](https://github.com/jump-dev/JuMP.jl) decision variables.
The commutativity of `T` with `*` is not assumed, even if it is the coefficient of a monomial of commutative variables.
However, commutativity of `T` and of the variables `+` is always assumed.
This allows to keep the terms sorted (Graded Lexicographic order is used) in polynomial and measure which enables more efficient operations.

Below is a simple usage example

```julia
julia> using DynamicPolynomials

julia> @polyvar x y # assigns x (resp. y) to a variable of name x (resp. y)
(x, y)

julia> p = 2x + 3.0x*y^2 + y # define a polynomial in variables x and y
y + 2.0x + 3.0xy²

julia> differentiate(p, x) # compute the derivative of p with respect to x
2.0 + 3.0y²

julia> differentiate.(p, (x, y)) # compute the gradient of p
(2.0 + 3.0y², 1.0 + 6.0xy)

julia> p((x, y)=>(y, x)) # replace any x by y and y by x
2.0y + x + 3.0x²y

julia> subs(p, y=>x^2) # replace any occurence of y by x^2
2.0x + x² + 3.0x⁵

julia> p(x=>1, y=>2) # evaluate p at [1, 2]
16.0
```
Below is an example with `@polyvar x[1:n]`

```julia
julia> n = 3;

julia> @polyvar x[1:n] # assign x to a tuple of variables x1, x2, x3
(Variable{DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}, Graded{LexOrder}}[x₁, x₂, x₃],)

julia> p = sum(x .* x) # compute the sum of squares
x₃² + x₂² + x₁²

julia> subs(p, x[1]=>2, x[3]=>3) # make a partial substitution
13 + x₂²

julia> A = reshape(1:9, 3, 3);

julia> p(x => A * vec(x))  # corresponds to dot(A*x, A*x), need vec to convert the tuple to a vector
194x₃² + 244x₂x₃ + 77x₂² + 100x₁x₃ + 64x₁x₂ + 14x₁²
```

The terms of a polynomial are ordered in increasing monomial order. The default
ordering is the graded lex order but it can be modified using the
`monomial_order` keyword argument of the `@polyvar` macro.
We illustrate this below by borrowing the example p. 59 of "Ideals, Varieties and Algorithms"
of Cox, Little and O'Shea:
```julia
julia> p(x, y, z) = 4x*y^2*z + 4z^2 - 5x^3 + 7x^2*z^2
p (generic function with 1 method)

julia> @polyvar x y z monomial_order = LexOrder
(x, y, z)

julia> p(x, y, z) 
4z² + 4xy²z + 7x²z² - 5x³

julia> @polyvar x y z
(x, y, z)

julia> p(x, y, z)
4z² - 5x³ + 4xy²z + 7x²z²

julia> @polyvar x y z monomial_order = Graded{Reverse{InverseLexOrder}}
(x, y, z)

julia> p(x, y, z)
4z² - 5x³ + 7x²z² + 4xy²z
```

Note that, when doing substitution, it is required to give the `Variable` ordering that is meant.
Indeed, the ordering between the `Variable` is not alphabetical but rather by order of creation
which can be undeterministic with parallel computing.
Therefore, this order cannot be used for substitution, even as a default (see [here](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl/issues/3) for a discussion about this).

[build-img]: https://github.com/JuliaAlgebra/DynamicPolynomials.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaAlgebra/DynamicPolynomials.jl/actions?query=workflow%3ACI
[codecov-img]: http://codecov.io/github/JuliaAlgebra/DynamicPolynomials.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaAlgebra/DynamicPolynomials.jl?branch=master

[zenodo-url]: https://doi.org/10.5281/zenodo.1203245
[zenodo-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.1203245.svg
