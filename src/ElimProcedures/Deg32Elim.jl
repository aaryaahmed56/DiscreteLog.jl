__precompile__()

module Deg32Elim

################################################################################
#
#  Submodules
#
################################################################################

include("../EllCrv/EllCrv.jl")
include("../EllCrv/Finite.jl")
include("../EllCrvModel/EllCrvModel.jl")

################################################################################
#
#  Imports
#
################################################################################

import AbstractAlgebra

# Main Verification Procedure
function verify(prime::UInt64)

    if prime = 2
        
        # Fix a field of prime order.
        Z::AbstractAlgebra.FinField = AbstractAlgebra.GF(prime)

        # Fix homogeneous coordinates [x:y:z] for an elliptic curve, 
        # where x, y, z are in Z.
        x = x::AbstractAlgebra.FinFieldElem ∈ Z
        y = y::AbstractAlgebra.FinFieldElem ∈ Z
        z = x::AbstractAlgebra.FinFieldElem ∈ Z

        # Fix coefficients A, B in Z for the Weierstrass form of the
        # elliptic curve, i.e. in homogenous coordinates [x:y:z],
        # EllCrv := x^3 + Ax^2 + B - x*y.

        A = A::AbstractAlgebra.FinFieldElem ∈ Z
        B = B::AbstractAlgebra.FinFieldElem ∈ Z



        # Instantiate Polynomial Rings.
        R0, A, B, xQ, yQ, x1, y1, z1 = AbstractAlgebra.PolynomialRing(Z, 
            "A",
            "B",
            "xQ",
            "yQ",
            "x1",
            "y1",
            "z1")
        
            # R1, aq1, aq, a1, a0, r = AbstractAlgebra.PolynomialRing(R0,
        #    "aq1",
        #    "aq",
        #    "a1",
        #    "a0",
        #    "r")

        
        # Instantiate Elliptic Curve.
        EllCrvWeier::EllCrv = EllCrvWeier
        
        
    end
    
end

# Add two points P, Q with projective coordinates
# P := [P[1]:P[2]:P[3]], Q := [Q[1]:Q[2]:Q[3]] 
function add(P::EllCrvPt, Q::EllCrvPt, A)
    x1 = (P[1])/(P[3])
    y1 = (P[2])/(P[3])
    x2 = (Q[1])/(Q[3])
    y2 = (Q[2])/(Q[3])

    X = x1 + x2
    Y = y1 + y2

    x3 = X + A + Y/X + (y1^2 + y2^2)/(x1^2 + x2^2)
    y3 = y1 + x3 + (x1 + x3)*Y/X

    g = gcd(Denominator(x3), Denominator(y3))

    return [(Numerator(x3)*Denominator(y3)) div g, 
    (Numerator(x3)*Denominator(x3)) div g, 
    (Denominator(x3)*Denominator(y3)) div g]

end