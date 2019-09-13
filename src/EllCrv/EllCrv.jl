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
using Markdown

################################################################################
#
#  Exports
#
################################################################################

export EllCrv, EllCrvPt

export base_field, disc, EllipticCurve, infinity, is_inifinite, 
    is_short, is_on_curve, j_invariant, Psi_polynomial, is_rational, 
    psi_poly_field, short_weierstrass_model, +, ×

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

mutable struct AbstractVariety{V}
end

mutable struct EllCrv{T} <: AbstractVariety
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
        
        # Short (Weierstrass) form of the elliptic curve,
        # i.e. y^2 = x^3 + ax + b.
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

#function isin(a, F::AbstractAlgebra.FieldElem{T}) where 
#    {T <: AbstractAlgebra.FieldElem}
#  iszero(a) && return true
#  F.comp || (!isone(gcd(denominator(a), prime(L))) && return false)
#  F.comp && ppio(denominator(a), prime(L))[1] != denominator(a.data) && return false
#  return true
#end

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

################################################################################
#
#  Parent
#
################################################################################

function parent_type(::Type{AbstractAlgebra.FieldElem{T}}) where T
end

################################################################################
#
#  Field access/ Field division
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

function a_invars(E::EllCrv)
  if isdefined(E, :a_invars)
    return [ deepcopy(z) for z in E.a_invars ]
  else
    t = (E.coeff[1], E.coeff[2], E.coeff[3], E.coeff[4], E.coeff[5])
    E.a_invars = t
    return t
  end
end

function b_invars(E::EllCrv)
  if isdefined(E, :long_b)
    return deepcopy(E.b_invars)
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
#   Basic Field Data, Quotients in Fields
#
###############################################################################

data(a::AbstractAlgebra.FieldElem) = a.data

numerator(a::AbstractAlgebra.FieldElem{T}) where 
{T <: AbstractAlgebra.FieldElement} = numerator(data(a))

denominator(a::AbstractAlgebra.FieldElem{T}) where 
{T <: AbstractAlgebra.FieldElem} = denominator(data(a))

prime(L::AbstractAlgebra.Field) = L.prime


################################################################################
#
#  Point at infinity
#
################################################################################

@doc Markdown.doc"""
    infinity(E::EllCrv) -> EllCrvPt

Creates the point at infinity.
"""
function infinity(E::EllCrv{T}) where T
  infi = EllCrvPt{T}(E)
  return infi
end


################################################################################
#
#  Test for inclusion
#
################################################################################

@doc Markdown.doc"""
    is_on_curve(E::EllCrv{T}, coords::Array{T, 1}) -> Bool

Returns true if `coords` defines a point  on E and false otherwise. The array
`coords` must have length 2.
"""
function is_on_curve(E::EllCrv{T}, coords::Array{T, 1}) where T
  length(coords) != 2 && error("Array must be of length 2")

  x = coords[1]
  y = coords[2]

  if E.short == true
    if y^2 == x^3 + (E.coeff[1])*x + (E.coeff[2])
      return true
    else
      return false
    end
  else
    if (y^2 + (E.coeff[1])*x*y + (E.coeff[3])*y ==
            x^3 + (E.coeff[2])*x^2+(E.coeff[4])*x + (E.coeff[5]))
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

@doc Markdown.doc"""
    disc(E::EllCrv{T}) -> T

Computes the discriminant of $E$.
"""
function disc(E::EllCrv{T}) where T
  if isdefined(E, :disc)
    return E.disc
  end

  if E.short == true
    # fall back to the formula for the long form
    # this should be done in a more clever way
    R = base_field(E)
    F = EllipticCurve([R(0), R(0), R(0), E.coeff[1], E.coeff[2]])
    d = disc(F)
    E.disc = d
    return d::T
  else
    a1 = E.coeff[1]
    a2 = E.coeff[2]
    a3 = E.coeff[3]
    a4 = E.coeff[4]
    a6 = E.coeff[5]

    (b2, b4, b6, b8) = b_invars(E)

    d = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    E.disc = d
    return d::T
  end
end

################################################################################
#
#  j-invariant
#
################################################################################

# p. 46 Washington, p. 72 Cohen
@doc Markdown.doc"""
    j(E::EllCrv{T}) -> T

Computes the j-invariant of $E$.
"""
function j_invariant(E::EllCrv{T}) where T
  if isdefined(E, :j)
    return E.j
  end

  if E.short == true

    R = base_field(E)
    F = EllipticCurve([R(0), R(0), R(0), E.coeff[1], E.coeff[2]])
    j = j_invariant(F)
    E.j = j
    return j::T
  else
    a1 = E.coeff[1]
    a2 = E.coeff[2]
    a3 = E.coeff[3]
    a4 = E.coeff[4]
    a6 = E.coeff[5]

    (b2, b4, b6, b8) = b_invars(E)
    #b2 = a1^2 + 4*a2
    #b4 = a1*a3 + 2*a4
    #b6 = a3^2 + 4*a6
    #b8 = a1^2*a6 - a1*a3*a4 + 4*a2*a6 + a2*a3^2 - a4^2

    c4 = b2^2 - 24*b4

    j = divexact(c4^3, disc(E))
    E.j = j
    return j::T
  end
end

################################################################################
#
#  Addition of Points
#
################################################################################

# washington p. 14, cohen p. 270
@doc Markdown.doc"""
    +(P::EllCrvPt, Q::EllCrvPt, coords::Array) -> EllCrvPt
Adds two points on an elliptic curve, whether affine 
or projective, or in the former case, whether in short 
Weierstrass form or long form.

** Not implemented in characteristic 2.
"""

# Add two points P, Q with projective coordinates
# P := [P[1]:P[2]:P[3]], Q := [Q[1]:Q[2]:Q[3]] or 
# affine coordinates P:= (P[1], P[2]), Q := (Q[1], Q[2]) 
function +(P::EllCrvPt{T}, Q::EllCrvPt{T}, coords::Array{T, 1}) where T
    parent(P) != parent(Q) && error("Points must be on the same curve.")

    # characteristic(base_field(parent(P))) == 2 &&
    #    error("Choose a characteristic for the field that is 
    #    at least >= 2.")

    if length(coords == 2)
        # Is either P or Q the point at infinity?
        if P.is_infinite
            return Q
        elseif Q.is_infinite
            return P
        end

        E = P.parent
        
        # Distinguish between long and short form with affine coordinates
        if E.short == true
            if P.affcoordx !=  Q.affcoordx
                m = divexact(Q.affcoordy - P.affcoordy, Q.affcoordx - P.affcoordx)
                x = m^2 - P.affcoordx - Q.affcoordx
                y = m*(P.affcoordx - x) - P.affcoordy
            elseif P.affcoordy != Q.affcoordy
                return infinity(E)
            elseif P.coordy != 0
                m = divexact(3*(P.affcoordx)^2 + (E.coeff[1]), 2*P.affcoordy)
                x = m^2 - 2*P.affcoordx
                y = m*(P.affcoordx - x) - P.affcoordy
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
            if P.affcoordx == Q.affcoordx
                if Q.affcoordy == -a1*P.affcoordx - a3 - P.affcoordy
                    return infinity(E)
                elseif P.affcoordy == Q.affcoordy
                    m = divexact(3*((P.affcoordx)^2) + 2*a2*P.affcoordx + a4 
                    - a1*P.coordy, 2*P.affcoordy + a1*P.affcoordx + a3)
                    x = -P.affcoordx - Q.affcoordx - a2 + a1*m
                    y = -P.affcoordy - m*(x - P.affcoordx) - a1*x - a3
                else
                    m = divexact(Q.affcoordy - P.affcoordy, Q.affcoordx - P.affcoordx)
                    x = -P.affcoordx - Q.affcoordx - a2 + a1*m + m^2
                    y = -P.affcoordy - m*(x - P.affcoordx) - a1*x - a3
                end
            else
                m = divexact(Q.affcoordy - P.affcoordy, Q.affcoordx - P.affcoordx)
                x = -P.affcoordx - Q.affcoordx - a2 + a1*m + m^2
                y = -P.affcoordy - m*(x - P.affcoordx) - a1*x - a3
            end

            Erg = E([x, y], false)

        end
        return Erg
    end
end

#    x1 = (P.projcoordx)/(P.projcoordz)
#    y1 = (P.projcoordy)/(P.projcoordz)
#    x2 = (Q.projcoordx)/(Q.projcoordz)
#    y2 = (Q.projcoordy)/(Q.projcoordz)
#
#    X = x1 + x2
#    Y = y1 + y2
#
#    x3 = X + A + Y/X + (y1^2 + y2^2)/(x1^2 + x2^2)
#    y3 = y1 + x3 + (x1 + x3)*Y/X
#
#    g = gcd(denominator(x3), denominator(y3))
#
#    return [(numerator(x3)*denominator(y3)) div g, 
#    (numerator(x3)*denominator(x3)) div g, 
#    (denominator(x3)*denominator(y3)) div g]
#
#end