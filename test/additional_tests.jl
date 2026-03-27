using Test
using DynamicPolynomials
using MultivariatePolynomials
import MultivariatePolynomials as MP

@testset "Division (divides)" begin
    @testset "Commutative divides" begin
        @polyvar x y z

        # Monomial divides itself
        @test divides(x^2, x^2)
        @test divides(x * y, x * y)

        # Lower power divides higher power
        @test divides(x, x^3)
        @test divides(x^2, x^5)

        # Higher power does not divide lower power
        @test !divides(x^3, x)
        @test !divides(x^2, x)

        # Constant (empty monomial) divides everything
        @test divides(x^0, x^2 * y^3)
        @test divides(x^0, y)

        # Multi-variable divides
        @test divides(x * y, x^2 * y^3)
        @test !divides(x^2 * y^3, x * y)

        # Variable not in divisor
        @test divides(x, x * y^2)
        @test !divides(x * y * z, x * y)

        # Disjoint variables
        @test !divides(z, x * y)
        @test divides(x, x * y * z)
    end

    @testset "Non-commutative divides error" begin
        @ncpolyvar a b
        @test_throws ErrorException divides(a, a * b)
    end
end

@testset "Complex variables" begin
    @testset "@complex_polyvar macro" begin
        @complex_polyvar z w
        @test z isa Variable
        @test w isa Variable
        @test !isreal(z)
        @test !isreal(w)
    end

    @testset "conj, real, imag of complex variables" begin
        @complex_polyvar z
        zc = conj(z)
        zr = real(z)
        zi = imag(z)

        # conj of complex variable
        @test isconj(zc)
        @test !isconj(z)

        # conj of conj returns original
        @test conj(zc) == z

        # real part properties
        @test isrealpart(zr)
        @test isreal(zr)

        # imag part properties
        @test isimagpart(zi)

        # ordinary_variable
        @test ordinary_variable(z) == z
        @test ordinary_variable(zc) == z
        @test ordinary_variable(zr) == z
        @test ordinary_variable(zi) == z

        # imag of conj returns -1 * imag_part
        zi_conj = imag(zc)
        @test coefficient(zi_conj) == -1
    end

    @testset "conj of real variable is identity" begin
        @polyvar x
        @test conj(x) == x
        @test isreal(x)
        @test real(x) == x
    end

    @testset "Complex variable in polynomials" begin
        @complex_polyvar z
        p = z^2 + conj(z)
        @test !isreal(p)

        # Evaluation
        val = 1 + 2im
        result = p(z => val)
        @test result == val^2 + conj(val)
    end
end

@testset "Differentiation" begin
    @polyvar x y

    @testset "Basic differentiation" begin
        # d/dx(x^3) = 3x^2
        @test differentiate(x^3, x) == 3x^2

        # d/dy(x^3) = 0
        @test differentiate(x^3, y) == 0

        # d/dx(x*y^2) = y^2
        @test differentiate(x * y^2, x) == y^2

        # d/dy(x*y^2) = 2*x*y
        @test differentiate(x * y^2, y) == 2 * x * y
    end

    @testset "Polynomial differentiation" begin
        p = 3x^3 + 2x * y - y^2 + 5
        @test differentiate(p, x) == 9x^2 + 2y
        @test differentiate(p, y) == 2x - 2y
    end

    @testset "Differentiation of constant polynomial" begin
        p = 5 * x^0
        # Differentiating a constant gives zero
        @test iszero(differentiate(p, x))
    end

    @testset "Higher-order differentiation" begin
        p = x^4 + 2x^3 * y
        dp = differentiate(p, x)
        ddp = differentiate(dp, x)
        @test ddp == 12x^2 + 12x * y
    end
end

@testset "Antidifferentiation edge cases" begin
    @polyvar x y

    @testset "Antidifferentiate constant polynomial" begin
        p = 3 * one(x)
        pi_x = antidifferentiate(p, x)
        @test differentiate(pi_x, x) == p
    end

    @testset "Antidifferentiate with respect to missing variable" begin
        p = x^2
        pi_y = antidifferentiate(p, y)
        @test pi_y == x^2 * y
    end

    @testset "Antidifferentiate and differentiate roundtrip" begin
        p = 2x^2 * y + 3y^3
        pi_y = antidifferentiate(p, y)
        @test differentiate(pi_y, y) == p
    end
end

@testset "Substitution" begin
    @polyvar x y z

    @testset "Partial substitution (subs)" begin
        p = x^2 + y + z
        q = subs(p, x => 2)
        @test q == 4 + y + z

        q2 = subs(p, y => 0, z => 0)
        @test q2 == x^2
    end

    @testset "Full evaluation" begin
        p = x^2 + y * z
        @test p(x => 1, y => 2, z => 3) == 7
    end

    @testset "Substitution with polynomial values" begin
        p = x^2 + y
        q = subs(p, x => y)
        @test q == y^2 + y
    end

    @testset "Eval error for missing variable" begin
        p = x + y
        err = ArgumentError(
            "Variable `y` was not assigned a value. Use `subs` to substitute only a subset of the variables.",
        )
        @test_throws err p(x => 1)
    end

    @testset "Evaluation of zero polynomial" begin
        p = x - x
        @test p(x => 5) == 0
    end

    @testset "Monomial evaluation" begin
        m = x^2 * y^3
        @test m(x => 2, y => 3) == 4 * 27
    end

    @testset "Variable evaluation" begin
        @test x(x => 42) == 42
    end

    @testset "Term evaluation" begin
        t = 3x^2
        @test t(x => 4) == 48
    end
end

@testset "Hash consistency" begin
    @polyvar x y

    @testset "Variable hash" begin
        @polyvar a
        b = a
        @test hash(a) == hash(b)

        # Different variables have different hashes (with high probability)
        @polyvar c
        @test hash(a) != hash(c)
    end

    @testset "Monomial hash" begin
        # Equal monomials have equal hashes
        m1 = x^2 * y
        m2 = x^2 * y
        @test hash(m1) == hash(m2)

        # Monomial equal to variable has same hash
        @test hash(x) == hash(Monomial(x))
        @test hash(x) == hash(DynamicPolynomials.Monomial([x, y], [1, 0]))

        # Constant monomial
        @test hash(x^0) == hash(1)
    end

    @testset "MonomialVector hash" begin
        mv1 = MonomialVector([x, y], [[1, 0], [0, 1]])
        mv2 = MonomialVector([x, y], [[1, 0], [0, 1]])
        @test hash(mv1) == hash(mv2)

        # Single element MonomialVector hashes like the monomial
        mv_single = MonomialVector([x], [[2]])
        @test hash(mv_single) == hash(x^2)
    end
end

@testset "MonomialVector operations" begin
    @polyvar x y z

    @testset "Degree functions" begin
        mv = MonomialVector([x, y], [[1, 0], [2, 1], [0, 3]])
        @test mindegree(mv) == 1
        @test maxdegree(mv) == 3
        mind, maxd = extdegree(mv)
        @test mind == 1
        @test maxd == 3
    end

    @testset "Empty MonomialVector degrees" begin
        mv = monomials([x, y], Int[])
        @test mindegree(mv) == 0
        @test maxdegree(mv) == 0
        @test extdegree(mv) == (0, 0)
    end

    @testset "Canonical form" begin
        # MonomialVector with unused variable
        mv = MonomialVector([x, y], [[1, 0], [2, 0]])
        c = DynamicPolynomials.canonical(mv)
        @test variables(c) == [x]
        @test length(c) == 2
    end

    @testset "Iteration" begin
        mv = MonomialVector([x, y], [[1, 0], [0, 1]])
        collected = collect(mv)
        @test length(collected) == 2
        @test collected[1] == x
        @test collected[2] == y
    end

    @testset "deleteat! and pop!" begin
        mv = MonomialVector([x, y], [[1, 0], [0, 1], [1, 1]])
        deleteat!(mv, 2)
        @test length(mv) == 2
        pop!(mv)
        @test length(mv) == 1
        @test mv[1] == x
    end

    @testset "isempty" begin
        mv = MonomialVector([x, y], Vector{Int}[])
        @test isempty(mv)
        @test length(mv) == 0
    end

    @testset "Copy" begin
        mv = MonomialVector([x, y], [[1, 0], [0, 1]])
        mv2 = copy(mv)
        @test mv == mv2
        # Verify it's a deep copy
        push!(mv.Z, [1, 1])
        @test length(mv) == 3
        @test length(mv2) == 2
    end

    @testset "Boolean indexing" begin
        mv = MonomialVector([x, y], [[1, 0], [0, 1], [1, 1]])
        sub = mv[[true, false, true]]
        @test length(sub) == 2
        @test sub[1] == x
        @test sub[2] == x * y
    end

    @testset "Integer indexing" begin
        mv = MonomialVector([x, y], [[1, 0], [0, 1], [1, 1]])
        sub = mv[[1, 3]]
        @test length(sub) == 2
    end
end

@testset "Polynomial operations" begin
    @polyvar x y

    @testset "Leading term/monomial/coefficient" begin
        p = 3x^2 + 2x * y + y^2
        @test leading_coefficient(p) == 1
        @test leading_monomial(p) == x^2 * y^0  # depends on ordering
        @test !iszero(leading_term(p))
    end

    @testset "Leading operations on zero polynomial" begin
        p = x - x
        @test leading_coefficient(p) == 0
        @test iszero(leading_term(p))
    end

    @testset "remove_leading_term" begin
        p = x^2 + 2x + 3
        q = remove_leading_term(p)
        @test nterms(q) == nterms(p) - 1
    end

    @testset "map_coefficients" begin
        p = 2x^2 + 3y
        q = map_coefficients(c -> 2c, p)
        @test q == 4x^2 + 6y
        # Original unchanged
        @test p == 2x^2 + 3y
    end

    @testset "map_coefficients zeroing out" begin
        p = 2x^2 + 3y + 4
        q = map_coefficients(c -> c > 3 ? 0 : c, p)
        @test q == 2x^2 + 3y
    end

    @testset "Polynomial negation" begin
        p = x^2 + 2y
        @test -p == -x^2 - 2y
        @test p + (-p) == 0
    end

    @testset "Polynomial copy" begin
        p = x^2 + y
        q = copy(p)
        @test p == q
    end

    @testset "remove_monomials" begin
        p = x^2 + 2x * y + 3y^2
        q = remove_monomials(p, [x^2, y^2])
        @test q == 2x * y
    end

    @testset "isapprox" begin
        p = 1.0x^2 + 2.0y
        q = 1.0x^2 + 2.0y + 1e-15
        @test p ≈ q
        @test !(p ≈ (1.0x^2 + 3.0y))
    end
end

@testset "Promote variables" begin
    @polyvar x y z

    @testset "Promote monomials with different variables" begin
        m1 = x^2
        m2 = y^3
        pm1, pm2 = promote_variables(m1, m2)
        @test variables(pm1) == variables(pm2)
        @test pm1 == m1
        @test pm2 == m2
    end

    @testset "Promote monomials with same variables" begin
        m1 = x^2 * y
        m2 = x * y^3
        pm1, pm2 = promote_variables(m1, m2)
        @test pm1 === m1  # no promotion needed
        @test pm2 === m2
    end
end

@testset "Monomial operations" begin
    @polyvar x y

    @testset "Monomial power" begin
        @test x^0 == 1
        @test x^1 == x
        @test x^3 == x * x * x
    end

    @testset "Monomial canonical" begin
        m = DynamicPolynomials.Monomial([x, y], [2, 0])
        c = DynamicPolynomials.canonical(m)
        @test variables(c) == [x]
        @test exponents(c) == [2]
        @test m == c
    end

    @testset "Monomial copy" begin
        m = x^2 * y^3
        m2 = copy(m)
        @test m == m2
        m2.z[1] = 0
        @test m != m2  # independent copy
    end

    @testset "Monomial conversion from 1" begin
        @polyvar x
        MT = monomial_type(x)
        m = convert(MT, 1)
        @test m == 1
        @test iszero(nvariables(m))
    end

    @testset "Monomial conversion error" begin
        @polyvar x
        MT = monomial_type(x)
        @test_throws ErrorException convert(MT, 2)
    end

    @testset "Monomial from variable" begin
        m = Monomial(x)
        @test m == x
        @test variables(m) == [x]
        @test exponents(m) == [1]
    end
end

@testset "Variable properties" begin
    @testset "iscomm" begin
        @polyvar x
        @test DynamicPolynomials.iscomm(typeof(x)) == true
        @ncpolyvar y
        @test DynamicPolynomials.iscomm(typeof(y)) == false
    end

    @testset "name and name_base_indices" begin
        @polyvar x
        @test name(x) == "x"
        base, indices = MP.name_base_indices(x)
        @test base == "x"
        @test indices == Int[]

        @polyvar a[1:3]
        base, indices = MP.name_base_indices(a[1])
        @test base == "a"
        @test indices == [1]
    end

    @testset "Variable ordering" begin
        @polyvar a b c
        # Variables created later have lower order (sorted reverse by creation)
        @test a > b
        @test b > c
    end

    @testset "Variable equality and identity" begin
        @polyvar x
        y = x
        @test x == y
        @polyvar x  # new x with different creation order
        @test y != x
    end
end

@testset "Polynomial from matrix (quadratic form)" begin
    @polyvar x y
    mv = monomials([x, y], 1)
    Q = [1 2; 3 4]
    p = polynomial(Q, mv)
    # p should be x^2 * Q[1,1] + (Q[1,2]+Q[2,1])*x*y + Q[2,2]*y^2
    # = x^2 + 5*x*y + 4*y^2
    @test p == x^2 + 5x * y + 4y^2
end

@testset "Monomial orderings in construction" begin
    @testset "LexOrder" begin
        @polyvar x y monomial_order = LexOrder
        monos = monomials([x, y], 2)
        @test length(monos) == 3
        @test issorted(monos)
    end

    @testset "Graded{LexOrder} (default)" begin
        @polyvar x y
        monos = monomials([x, y], 0:2)
        @test issorted(monos)
        @test first(monos) == x^0
    end

    @testset "InverseLexOrder" begin
        @polyvar x y monomial_order = InverseLexOrder
        monos = monomials([x, y], 2)
        @test length(monos) == 3
        @test issorted(monos)
    end
end

@testset "Arithmetic edge cases" begin
    @polyvar x y

    @testset "Addition of equal monomials" begin
        @test x + x == 2x
        @test x^2 + x^2 == 2x^2
    end

    @testset "Subtraction to zero" begin
        @test x - x == 0
        @test iszero(x^2 * y - x^2 * y)
    end

    @testset "Multiplication by zero" begin
        p = x^2 + y
        @test 0 * p == 0
        @test iszero(0 * p)
    end

    @testset "Multiplication by one" begin
        p = x^2 + y
        @test 1 * p == p
    end

    @testset "Mixed type arithmetic" begin
        p = x^2 + y
        q = 0.5 * p
        @test q == 0.5x^2 + 0.5y
        @test coefficient_type(q) == Float64
    end

    @testset "Polynomial + constant" begin
        p = x^2 + y
        @test p + 1 == x^2 + y + 1
        @test 1 + p == x^2 + y + 1
    end

    @testset "Polynomial * Polynomial" begin
        p = x + y
        q = x - y
        @test p * q == x^2 - y^2
    end
end

@testset "Non-commutative operations" begin
    @ncpolyvar x y

    @testset "Non-commutative multiplication order matters" begin
        @test x * y != y * x
    end

    @testset "Non-commutative polynomial arithmetic" begin
        p = x * y + y * x
        @test p != 2 * x * y
    end

    @testset "Non-commutative substitution" begin
        p = x * y * x
        q = subs(p, x => y)
        @test q == y^3
    end

    @testset "Non-commutative power" begin
        # x^2 for non-commutative is x*x
        @test x^2 == x * x
        @test (x * y)^1 == x * y
    end

    @testset "Non-commutative monomial construction" begin
        m = DynamicPolynomials.Monomial([x, y, x], [1, 2, 1])
        @test variables(m) == [x, y, x]
        @test exponents(m) == [1, 2, 1]
    end

    # Tests for PR #189: Fix multiplication of noncommutative variables
    # The fix ensures trailing zero-exponent variables are stripped from
    # the result of Monomial * Monomial multiplication.
    @testset "NC Monomial * Monomial trailing zeros (PR #189)" begin
        @ncpolyvar a b c

        # Multiplying monomials that have trailing zero exponents should
        # not carry those zeros into the result.
        m1 = DynamicPolynomials.Monomial([a, b], [1, 0])  # effectively just `a`, with trailing zero for `b`
        m2 = DynamicPolynomials.Monomial([a, b], [0, 1])  # effectively just `b`, with leading zero for `a`
        result = m1 * m2
        # Should be a*b without extra zero-exponent entries
        @test result == a * b
        @test all(!iszero, exponents(result))

        # Same variable at boundary: trailing zeros stripped when merging
        m3 = DynamicPolynomials.Monomial([a, b, c], [2, 0, 0])  # a^2 with trailing zeros
        m4 = DynamicPolynomials.Monomial([a, b, c], [0, 0, 3])  # c^3 with leading zeros
        result2 = m3 * m4
        @test result2 == a^2 * c^3
        @test all(!iszero, exponents(result2))

        # When last var of x matches first var of y (merge case)
        m5 = DynamicPolynomials.Monomial([a, b, c], [1, 1, 0])  # a*b with trailing zero
        m6 = DynamicPolynomials.Monomial([a, b, c], [0, 1, 1])  # b*c with leading zero
        result3 = m5 * m6
        @test result3 == a * b^2 * c
        @test all(!iszero, exponents(result3))

        # Multiplying by zero monomial (all zeros) returns the other
        m_zero = DynamicPolynomials.Monomial([a, b], [0, 0])
        m_val = DynamicPolynomials.Monomial([a, b], [1, 1])
        @test m_zero * m_val == a * b
        @test m_val * m_zero == a * b

        # Both sides have trailing zeros, different vars at boundary
        m7 = DynamicPolynomials.Monomial([a, b, c], [1, 0, 0])
        m8 = DynamicPolynomials.Monomial([a, b, c], [0, 0, 2])
        result4 = m7 * m8
        @test result4 == a * c^2
        @test all(!iszero, exponents(result4))

        # Verify the result works correctly in polynomial arithmetic
        p = result4 + a * b
        @test length(terms(p)) == 2
    end

    @testset "NC Monomial * Monomial with matching boundary variable" begin
        @ncpolyvar a b c
        # When the last nonzero var of x equals the first nonzero var of y,
        # exponents should be summed and trailing zeros dropped
        m1 = DynamicPolynomials.Monomial([a, b, c], [0, 2, 0])  # b^2 with zeros on both sides
        m2 = DynamicPolynomials.Monomial([a, b, c], [0, 3, 0])  # b^3 with zeros on both sides
        result = m1 * m2
        @test result == b^5
        @test variables(result) == [b]
        @test exponents(result) == [5]
    end

    @testset "NC Monomial * Monomial complex patterns" begin
        @ncpolyvar a b c
        # Multi-variable with internal structure and trailing zeros
        m1 = DynamicPolynomials.Monomial([a, b, a, b, c], [1, 2, 1, 0, 0])  # a*b^2*a with trailing zeros
        m2 = DynamicPolynomials.Monomial([a, b, c], [0, 0, 1])  # c with leading zeros
        result = m1 * m2
        @test all(!iszero, exponents(result))

        # Ensure result multiplied further still works
        result2 = result * a
        @test !iszero(result2)
    end
end

@testset "mergevars" begin
    @polyvar x y z

    @testset "Merge disjoint variable sets" begin
        vars, maps = DynamicPolynomials.mergevars([[x], [y]])
        @test length(vars) == 2
        @test x in vars
        @test y in vars
    end

    @testset "Merge overlapping variable sets" begin
        vars, maps = DynamicPolynomials.mergevars([[x, y], [y, z]])
        @test length(vars) == 3
        @test vars == [x, y, z]
    end

    @testset "Merge identical variable sets" begin
        vars, maps = DynamicPolynomials.mergevars([[x, y], [x, y]])
        @test vars == [x, y]
        @test maps[1] == maps[2]
    end
end

@testset "Monomial conj (commutative)" begin
    @complex_polyvar z
    m = z^2 * conj(z)^3
    mc = conj(m)
    # conj swaps z and conj(z)
    @test mc == conj(z)^2 * z^3
end

@testset "Zero polynomial properties" begin
    @polyvar x y
    p = zero(x + y)
    @test iszero(p)
    @test nterms(p) == 0
    @test variables(p) == [x, y]

    p_typed = zero(typeof(x + y))
    @test iszero(p_typed)
end

@testset "One polynomial properties" begin
    @polyvar x y
    p = one(x + y)
    @test !iszero(p)
    @test nterms(p) == 1
    @test p == 1
end
