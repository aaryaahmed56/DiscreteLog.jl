import AbstractAlgebra

################################################################################
#
#  Factorization Methods
#
################################################################################

# Square-free Factorization
function factor_sff(f::AbstractAlgebra.Generic.Poly{AbstractAlgebra.gfelem{Int32}})
    
    # Is monic?

    R = 1
    C = AbstractAlgebra.gcd(f, f')
    w = div(f, C)
end

# Distinct Degree Factorization
function factor_ddf(f::AbstractAlgebra.Generic.Poly{AbstractAlgebra.gfelem{Int32}})
    
    # Is square-free?
    
    i = 1
    k = AbstractAlgebra.base_ring(f)

    S, x = AbstractAlgebra.PolynomialRing(k, "x", cached = false)

    fstar = f

    while AbstractAlgebra.degree(fstar) >= 2*i
        r = x^(AbstractAlgebra.characteristic(k))^i - x
        g = AbstractAlgebra.gcd(fstar, r)
        if g != 1
            S = push!(S, (g, i))
            fstar = AbstractAlgebra.divexact(fstar, g)
        end
    end
    
    if fstar != 1
        S = push!(S, (fstar, AbstractAlgebra.degree(fstar)))
    elseif S == []
        return [(f, 1)]
    else
        return S
    end
end