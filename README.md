# Dynamic Polynomials

| **Build Status** | **References to cite** |
|:----------------:|:----------------------:|
| [![Build Status][build-img]][build-url] | [![DOI][zenodo-img]][zenodo-url] |
| [![Codecov branch][codecov-img]][codecov-url] | |

Sparse dynamic representation of multivariate polynomials that can be used with [MultivariatePolynomials](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl) (see the documentation there for more information).
Both commutative and non-commutative variables are supported.
The following types are defined:

* `Variable{C}`: A variable which is commutative with `*` when `C` is `true`. Commutative variables are created using the `@polyvar` macro, e.g. `@polyvar x y`, `@polyvar x[1:8]` and non-commutative variables are created likewise using the `@ncpolyvar` macro.
* `Monomial{C}`: A product of variables: e.g. `x*y^2`.
* `Term{C, T}`: A product between an element of type `T` and a `Monomial{C}`, e.g `2x`, `3.0x*y^2`.
* `Polynomial{C, T}`: A sum of `Term{C, T}`, e.g. `2x + 3.0x*y^2 + y`.

All common algebraic operations between those types are designed to be as efficient as possible without doing any assumption on `T`.
Typically, one imagine `T` to be a subtype of `Number` but it can be anything.
This is useful for example in the package [PolyJuMP](https://github.com/JuliaOpt/PolyJuMP.jl) where `T` is often an affine expression of [JuMP](https://github.com/JuliaOpt/JuMP.jl) decision variables.
The commutativity of `T` with `*` is not assumed, even if it is the coefficient of a monomial of commutative variables.
However, commutativity of `T` and of the variables `+` is always assumed.
This allows to keep the terms sorted (Graded Lexicographic order is used) in polynomial and measure which enables more efficient operations.

Below is a simple usage example

```julia
julia> using DynamicPolynomials

julia> @polyvar x y # assigns x (resp. y) to a variable of name x (resp. y)
(x, y)

julia> p = 2x + 3.0x*y^2 + y # define a polynomial in variables x and y
3.0xy² + 2.0x + y

julia> differentiate(p, x) # compute the derivative of p with respect to x
3.0y² + 2.0

julia> differentiate.(p, (x, y)) # compute the gradient of p
(3.0y² + 2.0, 6.0xy + 1.0)

julia> p((x, y)=>(y, x)) # replace any x by y and y by x
3.0x²y + x + 2.0y

julia> subs(p, y=>x^2) # replace any occurence of y by x^2
3.0x⁵ + x² + 2.0x

julia> p(x=>1, y=>2) # evaluate p at [1, 2]
16.0
```
Below is an example with `@polyvar x[1:n]`

```julia
julia> n = 3;

julia> @polyvar x[1:n] # assign x to a tuple of variables x1, x2, x3
(Variable{true}[x₁, x₂, x₃],)

julia> p = sum(x .* x) # compute the sum of squares
x₁² + x₂² + x₃²

julia> subs(p, x[1]=>2, x[3]=>3) # make a partial substitution
x₂² + 13

julia> A = reshape(1:9, 3, 3);

julia> p(x => A * vec(x))  # corresponds to dot(A*x, A*x), need vec to convert the tuple to a vector
14x₁² + 64x₁x₂ + 100x₁x₃ + 77x₂² + 244x₂x₃ + 194x₃²
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
