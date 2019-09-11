################################################################################
#
#          EllCrv/EllCrv.jl : Elliptic curves over general fields
#
# This file is part of Hecke.
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

export EllCrv, EllCrvPt, EllCrvDivisor

export base_field, discriminant, EllipticCurve, infinity, is_inifinite, 
    is_short, is_on_curve, j_invariant, Psi_polynomial, is_rational, 
    psi_poly_field, short_weierstrass_model, +, ×

################################################################################
#
#  Types
#
################################################################################

mutable struct EllCrv{T}
    base_field::AbstractAlgebra.Field
    short::Bool
    coeff::Array{T, 1}
    a_invars::Tuple{T, T, T, T, T}
    b_invars::Tuple{T, T, T, T}
    long_c::Array{T, 1}
    discriminant::T
    j::T

    torsion_points#::Array{EllCrvPt{T}, 1}
    torsion_structure#::Tuple{Array{Int, 1}, Array{EllCrvPt{T}, 1}}

    function EllCrv{T}(coeffs::Array{T, 1}, check::Bool= true) where {T}
        if length(coeffs) == 2
            if check
                disc = 4*coeffs[1]^3 + 27*coeffs[2]^2
                if disc != 0
                    E = new{T}()
                    E.short = true

                    E.coeff = [ deepcopy(z) for z ∈ coeffs]
                    E.base_field = parent(coeffs[1])

                else
                    error("Discriminant is zero.")
                end
            else
                E = new{T}()
                E.short = true
                E.coeff = [ deepcopy(z) for z ∈ coeffs]
                E.base_field = parent(coeffs[1])
            end
        elseif length(coeffs) == 5
            if check
                a1 = coeffs[1]
                a2 = coeffs[2]
                a3 = coeffs[3]
                a4 = coeffs[4]
                a6 = coeffs[5]

                b2 = a1^2 + 4*a2
                b4 = a1*a3 + 2*a4
                b6 = a3*2 + 4*a6
                b8 = a1^2*a6 - a1*a3*a4 + 4*a2*a6 + a2*a3*2 - a4*2

                disc = (-b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6)

                if disc != 0
                    E = new{T}()
                    E.coeff = [ deepcopy(z) for z ∈ coeffs]
                    E.short = false
                    E.b_invars = (b2, b4, b6, b8)
                    E.a_invars = (a1, a2, a3, a4, a6)
                    E.discriminant = disc
                    E.base_field = parent(coeffs[1])
                else
                    error("Discriminant is zero.")
                end
            else
                E = new{T}()
                E.short = false
                E.coeff = [ deepcopy(z) for z ∈ coeffs]
                E.base_field = parent(coeffs[1])
            end
        else
            error("Length of coefficient array must be 2 or 5.")
        end
        return E 
    end
end

mutable struct EllCrvPt{T}
    # Accommodates affine as well as projective
    # coordinates for any point P on the elliptic curve.
    affcoordx::T
    affcoordy::T
    projcoordx::T 
    projcoordy::T
    projcoordz::T
    is_inifinite::Bool
    parent::EllCrv{T}

    function EllCrvPt{T}(E::EllCrv{T}, coords::Array{T, 1}, 
        check::Bool = true) where {T}
        if check
            if is_on_curve(E, coords)
                if length(coords) < 2
                    error("Point must have at least two coordinates.")
                elseif length(coords) == 2
                    # Point is presumed to have affine oords.
                    P = new{T}(coords[1], coords[2], false, E)
                    return P
                elseif length(coords) == 3
                    # Point is presumed to have projective coords.
                    P = new{T}(coords[1], coords[2], coords[3], false, E)
                    return P 
            else
                error("Point is not on curve.")
            end
        end
    end

    function EllCrvPt{T}(E::EllCrv{T}) where {T}
        z = new{T}()
        z.parent = E 
        z.is_inifinite = true
        return z
    end
end

mutable struct EllCrvDivisor{T, P_1,...,P_s}
    coeff::Array{T, 1}
    s = length(coeff)
    points::Tuple{P_1::EllCrvPt{T},...,P_s::EllCrvPt{T}}
    func_field_without_zero::AbstractAlgebra.Field
    rat_func::AbstractAlgebra.FieldElem ∈ func_field_without_zero
    degree::Int
    is_associated::Bool
    is_effective::Bool

    function EllCrvDivisor{T, P_1,...,P_s}(rat::AbstractAlgebra.FieldElem, 
        coeffs::Array{T, 1}, points::Tuple{P_1::EllCrvPt{T},...,P_s::EllCrvPt{T}}, 
        check::Bool = true) where {T, P_1,...,P_s}
        if check
            if is_rational(rat, coeffs)
                ECD = new{T, P_1,...,P_s}()
                ECD.is_associated = true
            end
        end
    end
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllipticCurve(x::Array{T, 1}, check::Bool = true) where T
    E = EllCrv{T}(x, check)
    return E
end

function(E::EllCrv{T})(coords::Array{S, 1}, check::Bool = true) where {S, T}
    if length(coords) < 2
        error("Need at least two coordinates.")
    elseif length(coords) == 2
        print("Point is affine.")
    elseif length(coords) == 3
        print("Point is homogeneous.")
    end

    if S == T
        parent(coords[1]) != base_field(E) &&
            error("The elliptic curve and the point must be defined over the same 
            field.")
            return EllCrvPt{T}(E, coords, check)
    else
        return EllCrvPt{T}(E, map(base_field(E), coords), check)
    end
end

#function EllipticCurveDiv(x::Array{T}, y::Tuple{P_1,...,P_s}, check::Bool = true) 
#    where {P_1,...,P_s, T}
#    EDiv = EllCrvDivisor{T, P_1,...,P_s, check}
#    return EDiv
#end

################################################################################
#
#  Field access
#
################################################################################

function base_field(E::EllCrv{T}) where T
    return E.base_field::parent_type(T)
end

function Base.deepcopy_internal(E::EllCrv, dict::IdDict)
    return EllipticCurve(E.coeff)
end

function parent(P::EllCrvPt)
    return P.parent
end

function is_finite(P::EllCrvPt)
    return !P.is_inifinite
end

function is_infinite(P::EllCrvPt)
    return P.is_infinite
end

function is_short(E::EllCrv)
    return E.short
end

function a


################################################################################
#
#  Addition of Points
#
################################################################################

# washington p. 14, cohen p. 270
@doc Markdown.doc"""
    +(P::EllCrvPt, Q::EllCrvPt) -> EllCrvPt
Adds two points on an elliptic curve.
does not work in characteristic 2
"""

# Add two points P, Q with projective coordinates
# P := [P[1]:P[2]:P[3]], Q := [Q[1]:Q[2]:Q[3]] 
function +(P::EllCrvPt{T}, Q::EllCrvPt{T}) where T
    parent(P) != parent(Q) && error("Points must be on the same curve.")

    characteristic(base_field(parent(P))) == 2 &&
        error("Choose a characteristic for the field that is 
        at least >= 2.")

      # Is P = infinity or Q = infinity?
    if P.is_infinite
        return Q
    elseif Q.is_infinite
        return P
    end

    E = P.parent

    x1 = (P.projcoordx)/(P.projcoordz)
    y1 = (P.projcoordy)/(P.projcoordz)
    x2 = (Q.projcoordx)/(Q.projcoordz)
    y2 = (Q.projcoordy)/(Q.projcoordz)

    X = x1 + x2
    Y = y1 + y2

    x3 = X + A + Y/X + (y1^2 + y2^2)/(x1^2 + x2^2)
    y3 = y1 + x3 + (x1 + x3)*Y/X

    g = gcd(Denominator(x3), Denominator(y3))

    return [(Numerator(x3)*Denominator(y3)) div g, 
    (Numerator(x3)*Denominator(x3)) div g, 
    (Denominator(x3)*Denominator(y3)) div g]

end

function +(P::EllCrvPt{T}, Q::EllCrvPt{T}) where T
  parent(P) != parent(Q) && error("Points must live on the same curve")

  characteristic(base_field(parent(P))) == 2 &&
      error("Not yet implemented in characteristic 2")



  E = P.parent

  # Distinguish between long and short form
  if E.short == true
    if P.coordx != Q.coordx
        m = divexact(Q.coordy - P.coordy, Q.coordx - P.coordx)
        x = m^2 - P.coordx - Q.coordx
        y = m * (P.coordx - x) - P.coordy
    elseif P.coordy != Q.coordy
        return infinity(E)
    elseif P.coordy != 0
        m = divexact(3*(P.coordx)^2 + (E.coeff[1]), 2*P.coordy)
        x = m^2 - 2*P.coordx
        y = m* (P.coordx - x) - P.coordy
    else
        return infinity(E)
    end

    Erg = E([x, y], false)

  else
    a1 = E.coeff[1]
    a2 = E.coeff[2]
    a3 = E.coeff[3]
    a4 = E.coeff[4]
    a6 = E.coeff[5]

    # Use [Cohen, p. 270]
    if P.coordx == Q.coordx
      if Q.coordy == -a1*P.coordx - a3 - P.coordy # then P = -Q
        return infinity(E)
      elseif P.coordy == Q.coordy # then P = Q
        m = divexact(3*((P.coordx)^2) + 2*a2*P.coordx + a4 - a1*P.coordy, 2*P.coordy + a1*P.coordx + a3)
        x = -P.coordx - Q.coordx - a2 + a1*m + m^2
        y = -P.coordy - m*(x - P.coordx) - a1*x - a3
      else # then P != +-Q
        m = divexact(Q.coordy - P.coordy, Q.coordx - P.coordx)
        x = -P.coordx - Q.coordx - a2 + a1*m + m^2
        y = -P.coordy - m*(x - P.coordx) - a1*x - a3
      end
    else # now P != +-Q
      m = divexact(Q.coordy - P.coordy, Q.coordx - P.coordx)
      x = -P.coordx - Q.coordx - a2 + a1*m + m^2
      y = -P.coordy - m*(x - P.coordx) - a1*x - a3
    end

    Erg = E([x, y], false)

  end
  return Erg
end