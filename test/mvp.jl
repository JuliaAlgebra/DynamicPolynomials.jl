using Pkg
import MultivariatePolynomials
const mvp_test =
    joinpath(dirname(pathof(MultivariatePolynomials)), "..", "test")
const Mod = DynamicPolynomials
const MP = MultivariatePolynomials
include(joinpath(mvp_test, "utils.jl"))
include(joinpath(mvp_test, "commutativetests.jl"))
include(joinpath(mvp_test, "noncommutativetests.jl"))
