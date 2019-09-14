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

export AbstractVariety, EllCrv, EllCrvPt
export base_field, disc, EllipticCurve, infinity, is_infinite, is_on_curve, 
j_invariant, +, ×

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

mutable struct AbstractVariety{V} end

mutable struct EllCrv{T} <: AbstractVariety
    base_field::AbstractAlgebra.Field
    coeff::Array{T, 1}
    a_invars::Tuple{T,...}
    b_invars::Tuple{T,...}
    long_c::Array{T, 1}
    j::T

    function EllCrv{T}(coeffs::Array{T, 1}, check::Bool = true) where T

        if check
            disc = 4*coeffs[1]^3 + 27*coeffs[2]^2
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

    function EllCrvPt{T}(E::EllCrv{T}) where T
        z = new{T}()
        z.parent = E
        z.is_infinite = true
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

@inline function parent_type(::Type{AbstractAlgebra.FieldElem{T}}) where T end
@inline function P.parent(P::EllCrvPt) return P.parent end
@inline function P.is_infinite(P::EllCrvPt) return P.is_infinite end

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
#   Basic Field Data, Quotients in Fields
#
###############################################################################

data(a::AbstractAlgebra.FieldElem) = a.data

numerator(a::AbstractAlgebra.FieldElem{T}) where 
{T <: AbstractAlgebra.FieldElem} = numerator(data(a))

denominator(a::AbstractAlgebra.FieldElem{T}) where 
{T <: AbstractAlgebra.FieldElem} = denominator(data(a))

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
            m = divexact(Q.coords[2] - P.coords[2], Q.coords[1] - P.coords[1])
            x = m^2 - P.coords[1] - Q.coords[1]
            y = m*(P.coords[1] - x) - P.coords[2]
        elseif P.coords[2] != Q.coords[2]
            return infinity(E)
        elseif P.coords[2] != 0
            m = divexact(3*(P.coords[1])^2 + (E.coeff[1]), 2*(P.coords[2]))
            x = m^2 - 2*(P.coords[1])
            y = m*(P.coords[1] - x) - P.coords[2]
        else
            return infinity(E)
        end

        Erg = E([x, y], false)
    end
    return Erg
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