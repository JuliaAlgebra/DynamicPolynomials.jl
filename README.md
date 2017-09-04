# Dynamic Polynomials

| **PackageEvaluator** | **Build Status** |
|:--------------------:|:----------------:|
| [![][pkg-0.6-img]][pkg-0.6-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][pkg-0.7-img]][pkg-0.7-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

Sparse dynamic representation of multivariate polynomials that can be used with [MultivariatePolynomials](https://github.com/blegat/MultivariatePolynomials.jl) (see the documentation there for more information).
Both commutative and non-commutative variables are supported.
The following types are defined:

* `PolyVar{C}`: A variable which is commutative with `*` when `C` is `true`. Commutative variables are created using the `@polyvar` macro, e.g. `@polyvar x y`, `@polyvar x[1:8]` and non-commutative variables are created likewise using the `@ncpolyvar` macro.
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
@polyvar x y # assigns x (resp. y) to a variable of name x (resp. y)
p = 2x + 3.0x*y^2 + y
@test differentiate(p, x) # compute the derivative of p with respect to x
@test differentiate.(p, (x, y)) # compute the gradient of p
@test p((x, y)=>(y, x)) # replace any x by y and y by x
@test subs(p, y=>x^2) # replace any occurence of y by x^2
@test p(x=>1, y=>2) # evaluate p at [1, 2]
```
Below is an example with `@polyvar x[1:n]`
```julia
n = 3
A = rand(n, n)
@polyvar x[1:n] # assign x to a tuple of variables x1, x2, x3
p = sum(x .* x) # x_1^2 + x_2^2 + x_3^2
subs(p, x[1]=>2, x[3]=>3) # x_2^2 + 13
p(x=>A*vec(x)) # corresponds to dot(A*x, A*x), need vec to convert the tuple to a vector
```
Note that, when doing substitution, it is required to give the `PolyVar` ordering that is meant.
Indeed, the ordering between the `PolyVar` is not alphabetical but rather by order of creation
which can be undeterministic with parallel computing.
Therefore, this order cannot be used for substitution, even as a default (see [here](https://github.com/blegat/MultivariatePolynomials.jl/issues/3) for a discussion about this).

[pkg-0.6-img]: http://pkg.julialang.org/badges/DynamicPolynomials_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=DynamicPolynomials
[pkg-0.7-img]: http://pkg.julialang.org/badges/DynamicPolynomials_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=DynamicPolynomials

[build-img]: https://travis-ci.org/blegat/DynamicPolynomials.jl.svg?branch=master
[build-url]: https://travis-ci.org/blegat/DynamicPolynomials.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/wu5dnoq4x3jvjft8?svg=true
[winbuild-url]: https://ci.appveyor.com/project/blegat/dynamicpolynomials-jl
[coveralls-img]: https://coveralls.io/repos/github/blegat/DynamicPolynomials.jl/badge.svg?branch=master
[coveralls-url]: https://coveralls.io/github/blegat/DynamicPolynomials.jl?branch=master
[codecov-img]: http://codecov.io/github/blegat/DynamicPolynomials.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/blegat/DynamicPolynomials.jl?branch=master
