const mvp_test = joinpath(Pkg.dir("MultivariatePolynomials"), "test")
const Mod = DynamicPolynomials
include(joinpath(mvp_test, "commutativetests.jl"))
