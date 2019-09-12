__precompile__()

module Deg32Elim

################################################################################
#
#  Submodules
#
################################################################################

include("../EllCrv/EllCrv.jl")
include("../EllCrv/Fields.jl")
include("../EllCrvModel/EllCrvModel.jl")

################################################################################
#
#  Imports
#
################################################################################

import AbstractAlgebra

using Markdown

@doc Markdown.doc"""
    verify(prime::UInt128) -> ...

An aside from the proof of Proposition 4.4, verifying 
the non-vanishing of a constructed determinant given 
a divisor D1 on E.
"""
function verify(prime::UInt128)

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

        coeff::Array{Z, 1} = {A, B}



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

        
        # Instantiate short Weierstrass Elliptic Curve.
        check = true
        EllCrvWeier = EllipticCurve(coeff, check)
        
        
    end
    
end