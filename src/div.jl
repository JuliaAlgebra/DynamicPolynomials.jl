# _div(a, b) assumes that b divides a
function _div(m1::Monomial{true}, m2::Monomial{true})
    w, updatez = multdivmono(m1.vars, m2, -)
    Monomial{true}(w, updatez(m1.z))
end
