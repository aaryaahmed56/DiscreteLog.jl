################################################################################
#
#          EllCrv/EllCrv.jl : Elliptic curves over general fields
#
# Copyright (c) 2015, 2016: Claus Fieker, Tommy Hofmann
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# (C) 2016 Tommy Hofmann
# (C) 2016 Robin Ammon
# (C) 2016 Sofia Brenner
#
################################################################################


################################################################################
#
#  Imports
#
################################################################################

import AbstractAlgebra

################################################################################
#
#  Exports
#
################################################################################

export AbstractVariety, EllCrv, EllCrvPt
export base_field, disc, EllipticCurve, infinity, generic_point, is_infinite, 
is_on_curve, j_invariant, +, *

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

mutable struct AbstractVariety{V} end

mutable struct AbstractCurve{C} end

mutable struct EllCrv{T} <: AbstractCurve where T <: AbstractAlgebra.FieldElem
    base_field::AbstractAlgebra.Field
    coeff::Array{T, 1}
    a_invars::Tuple{T,...}
    b_invars::Tuple{T,...}
    long_c::Array{T, 1}
    j::T

    function EllCrv{T}(coeffs::Array{T, 1}, check::Bool = true) where T

        if check
            disc = -16*(4*coeffs[1]^3 + 27*coeffs[2]^2)
            if disc != 0
                E = new{T}()
                E.coeff = [ deepcopy(z) for z ∈ coeffs]
                E.base_field = parent(coeffs[1])
            else
                error("Discriminant is zero.")
            end
        else
            E = new{T}()
            E.coeff = [ deepcopy(z) for z ∈ coeffs]
            E.base_field = parent(coeffs[1])
        end
        return E
    end
end

mutable struct EllCrvPt{T}
    coord::Array{T, 1}
    is_infinite::Bool
    is_generic::Bool
    parent::EllCrv{T}

    function EllCrvPt{T}(E::EllCrv{T}, coords::Array{T, 1}, 
        check::Bool = true) where T
        if check
            if is_on_curve(E, coords)
                if length(coords) < 2
                    error("Point must have at least two coordinates.")
                elseif length(coords) == 2
                    P = new{T}(coords[1], coords[2], false, E)
                elseif length(coords) == 3
                    P = new{T}(coords[1], coords[2], coords[3], false, E)
                end
                return P
            else
                error("Point is not on curve.")
            end
        end
    end

    # Point at infinity
    function EllCrvPt{T}(E::EllCrv{T}) where T
        z = new{T}()
        z.parent = E
        z.is_infinite = true
        return z
    end

    # Generic Point
    function EllCrvPt{T}(E::EllCrv{T}) where T
        z = new{T}()
        z.parent = E
        z.is_generic = true
        return z
    end
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllipticCurve{T}(coeffs::Array{T, 1}, check::Bool = true) where T
    E = EllCrv{T}(coeffs, check)
    return E
end

function(E::EllCrv{T})(coords::Array{S, 1}, check::Bool = true) where {S, T}
    if length(coords) < 2
        error("Point must have at least two coordinates.")
    elseif length(coords) == 2
        print("Point is affine.")
    elseif length(coords) == 3
        print("Point is projective.")
    end

    if S == T 
        parent(coords[1]) != base_field(E) &&
            error("The elliptic curve and point must be defined over the same
            field.")
            return EllCrvPt{T}(E, coords, check)
    else
        return EllCrvPt{T}(E, map(base_field(E), coords), check)
    end
end

################################################################################
#
#  Field access
#
################################################################################

@inline function parent_type(::Type{AbstractAlgebra.FieldElem}) end
@inline function P.parent(P::EllCrvPt) return P.parent end
@inline function P.is_infinite(P::EllCrvPt) return P.is_infinite end
@inline function P.is_generic(P::EllCrvPt) return P.is_generic end

function base_field(E::EllCrv{T}) where T
    return E.base_field::parent_type(T)
end

function a_invars(E::EllCrv)
    if isdefined(E, :a_invars)
        return [ deepcopy(z) for z ∈ E.a_invars ]
    else
        t = (E.coeff[1], E.coeff[2], E.coeff[3], E.coeff[4], E.coeff[5])
        E.a_invars = t
        return t
    end
end

function b_invars(E::EllCrv)
    if isdefined(E, :b_invars)
        return [ deepcopy(z) for z ∈ E.b_invars]
    else
        a1 = E.coeff[1]
        a2 = E.coeff[2]
        a3 = E.coeff[3]
        a4 = E.coeff[4]
        a6 = E.coeff[5]

        b2 = a1^2 + 4*a2
        b4 = a1*a3 + 2*a4
        b6 = a3^2 + 4*a6
        b8 = a1^2*a6 - a1*a3*a4 + 4*a2*a6 + a2*a3^2 - a4^2

        E.b_invars = (b2, b4, b6, b8)
        return (b2, b4, b6, b8)
    end
end

###############################################################################
#
#   Basic Field Data
#
###############################################################################

prime(L::AbstractAlgebra.Field) = L.prime

################################################################################
#
#  Point at infinity
#
################################################################################

function infinity(E::EllCrv{T}) where T 
    infi = EllCrvPt{T}(E)
    return infi
end

################################################################################
#
#  Generic Point
#
################################################################################

function generic_point(E::EllCrv{T}) where T 
    gen = EllCrvPt{T}(E)
    return generic_point
end

################################################################################
#
#  Test for inclusion
#
################################################################################

function is_on_curve(E::EllCrv{T}, coords::Array{T, 1}) where T 
    length(coords) != 2 && error("Point must have at least two coordinates.")

    if length(coords) == 2
        x = coords[1]
        y = coords[2]
        
        if y^2 == x^3 + (E.coeff[1])*x + (E.coeff[2])
            return true
        else
            return false
        end
    elseif length(coords) == 3
        x = coords[1]
        y = coords[2]
        z = coords[3]

        if y^2*z == x^3 + (E.coeff[1])*x*z^2 + (E.coeff[2])*z^3
            return true
        else
            return false
        end
    end
end

################################################################################
#
#  Discriminant
#
################################################################################

function disc(E::EllCrv{T}) where T 
    if isdefined(E, :disc)
        return E.disc
    end

    R = base_field(E)
    F = EllipticCurve([R(0), R(0), R(0), E.coeff[1], E.coeff[2]])
    d = disc(F)

    E.disc = d 
    return d::T 
end

################################################################################
#
#  j-invariant
#
################################################################################

function j_invariant(E::EllCrv{T}) where T
    if isdefined(E, :j)
        return E.j
    end

    R = base_field(E)
    F = EllipticCurve([R(0), R(0), R(0)], E.coeff[1], E.coeff[2])
    j = j_invariant(F)
    E.j = j
    return j::T
end

################################################################################
#
#  Scalar multiplication
#
################################################################################

function *(n::Int, P::EllCrvPt)
    C = P
    E = P.parent
    B = infinity(E)
    F = base_field(E)

    if AbstractAlgebra.characteristic(F) == 2
        # Implement Fast Scalar Mult for Binary Fields.
        # https://eprint.iacr.org/2017/840.pdf
    end 
    if n >= 0
        a = n
    else
        a = -n
    end
  
    while a != 0
      if mod(a,2) == 0
        a = div(a,2)
        C = C + C
      else
        a = a - 1
        B = B + C
      end
    end
  
    if n < 0
      B = -B
    end
  
    return B
  end

################################################################################
#
#  Addition of Points
#
################################################################################

# Add two points P, Q with projective coordinates
# P := [P[1]:P[2]:P[3]], Q := [Q[1]:Q[2]:Q[3]] or 
# affine coordinates P:= (P[1], P[2]), Q := (Q[1], Q[2])

function +(P::EllCrvPt, Q::EllCrvPt, coords::Array{T, 1}) where T 
    parent(P) != parent(Q) && error("Points must live on the same curve.")

    if length(coords) == 2
        
        # Is either P or Q the point at infinity?
        if P.is_infinite
            return Q
        elseif Q.is_infinite
            return P
        end

        E = P.parent

        if P.coords[1] != Q.coords[1]
            m = AbstractAlgebra.divexact(Q.coords[2] - P.coords[2], Q.coords[1] - P.coords[1])
            x = m^2 - P.coords[1] - Q.coords[1]
            y = m*(P.coords[1] - x) - P.coords[2]
        elseif P.coords[2] != Q.coords[2]
            return infinity(E)
        elseif P.coords[2] != 0
            m = AbstractAlgebra.divexact(3*(P.coords[1])^2 + (E.coeff[1]), 2*(P.coords[2]))
            x = m^2 - 2*(P.coords[1])
            y = m*(P.coords[1] - x) - P.coords[2]
        else
            return infinity(E)
        end

        Erg = E([x, y], false)
    end
    return Erg
end

################################################################################
#
#  Division Polynomials (for SEA)
#
################################################################################

# Make this dynamic
function division_polynomialE(E::EllCrv, n::Int, x = nothing, y = nothing)
    A = numerator(E.coeff[1])
    B = numerator(E.coeff[2])

    if x === nothing
        Z = AbstractAlgebra.ZZ
        Zx, _x = AbstractAlgebra.PolynomialRing(Z, "x")
        Zxy, y = AbstractAlgebra.PolynomialRing(Zx, "y")
    else
        Zxy = AbstractAlgebra.parent(x)
    end

    if n == 1
        return one(Zxy)
    elseif n == 2
        return 2*y
    elseif n == 3
        return 3*x^4 + 6*A*x^2 + 12*B*x - A^2
    elseif n == 4
        return 4*y*(x^6 + 5*A*x^4 + 20*B*x^3 - 5*A^2*x^2 - 4*A*B*x - 8*B^2 - A^3)
    elseif mod(n, 2) == 0 && n >= 3
        m = div(n, 2)
        return AbstractAlgebra.divexact(division_polynomialE(E, m, x, y), 2*y)
        *(division_polynomialE(E, m + 2, x, y)*division_polynomialE(E, m - 1, x, y)^2
        - division_polynomialE(E, m - 2, x, y)*division_polynomialE(E, m + 1, x, y)^2)
    else m = div(n - 1, 2)
        return division_polynomialE(E, m + 2, x, y)*division_polynomialE(E, m, x, y)^3 
        - division_polynomialE(E, m - 1, x, y)*division_polynomialE(E, m + 1, x, y)^3
    end
end