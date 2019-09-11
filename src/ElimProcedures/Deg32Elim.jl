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
        EllCrvWeier::EllCrv{T} = EllCrvWeier
        
        
    end
    
end