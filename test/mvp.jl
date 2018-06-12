using Compat.Pkg
const mvp_test = joinpath(Pkg.dir("MultivariatePolynomials"), "test")
const Mod = DynamicPolynomials
const MP = MultivariatePolynomials
include(joinpath(mvp_test, "utils.jl"))
include(joinpath(mvp_test, "commutativetests.jl"))
include(joinpath(mvp_test, "noncommutativetests.jl"))
