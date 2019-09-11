__precompile__()

module Deg32Elim

include("../EllCrv/EllCrv.jl")
include("../EllCrv/Finite.jl")
include("../EllCrvModel/EllCrvModel.jl")



using Pkg
using AbstractAlgebra


function verify(prime::UInt64)

    if prime = 2
        
        # Fix a field of prime order.
        Z::FinField = GF(prime)

        # Fix homogeneous coordinates [x:y] for an elliptic curve, 
        # where x, y are in Z.
        x = x::FinFieldElem ∈ Z
        y = x::FinFieldElem ∈ Z

        # Fix coefficients A, B in Z for the Weierstrass form of the
        # elliptic curve, i.e. in homogenous coordinates [x: y],
        # EllCrv := x^3 + Ax^2 + B - x*y.

        A = A::FinFieldElem ∈ Z
        B = B::FinFieldElem ∈ Z



        # Instantiate Polynomial Rings.
        R0, A, B, xQ, yQ, x1, y1, z1 = PolynomialRing(Z, 
            "A",
            "B",
            "xQ",
            "yQ",
            "x1",
            "y1",
            "z1")
        R = AbstractAlgebra.PolynomialRing(R0, 6)

        
        # Instantiate Elliptic Curve.
        EllCrvWeier::EllCrv{T} = EllCrv{T}()
    end
end



function add(P::AbstractVector, Q:: AbstractVector)

            x1 = (R!P[1])