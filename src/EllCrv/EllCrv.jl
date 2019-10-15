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

export AbstractVariety, AbstractCurve, EllCrv, EllCrvPt
export base_field, disc, EllipticCurve, infinity, isinfinite, isshort, 
issimp, ison_curve, j_invariant, ⊟, ⊞, ×

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

# Ancestral abstract type -- Abstract k-variety where k field 
# (Type::Union{AbstractAlgebra.GFField{Int64}, AbstractAlgebra.Field}).
abstract type AbstractVariety{T} end

abstract type AbstractCurve{T} <: AbstractVariety{T} end

mutable struct EllCrv{T} <: AbstractCurve{T}
    base_field::AbstractAlgebra.GFField{Int64}
    coeff::Array{T, 1}

    # General Weierstrass form invariants
    # for both projective and affine coordinates, i.e.
    # 
    # Affine (x, y):
    # y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    # 
    # Projective (X:Y:Z):
    # Y^2*Z + a1*X*Y*Z + a3*Y*Z^2 = X^3 + a2*X^2*Z + a4*X*Z^2 + a6*Z^3
    a_invars::Tuple{Vararg{T, N} where N}

    # Perform change of variables
    # phi: y --> 2*y + a1*x + a3
    # to determine b2, b4, b6.
    # b8 is for discriminant.
    b_invars::Tuple{Vararg{T, M} where M}
    c_invars::Array{T, 1}
    
    # Of form
    # Affine (x, y): y^2 = x^3 + A*x + B
    # Projective (X:Y:Z) = Y^2*Z = X^3*Z + A*X*Z^2 + B*Z^3
    # i.e. a1 = 0, a2 = 0, a3 = 0, a4 = A, a6 = B
    short::Bool
    
    # Of form
    # Affine (x, y): y^2 + x*y = x^3 + A*x^2 + B
    # Projective (X:Y:Z) = Y^2*Z + X*Y*Z = X^3 + A*X^2*Z + B*Z^3
    # i.e. a1 = 1, a2 = A, a3 = 0, a4 = 0, a6 = B
    simp::Bool
    disc::T
    j::T

    function EllCrv{T}(coeffs::Array{T, 1}, check::Bool = true) where T
        if check
            a1 = coeffs[1]
            a2 = coeffs[2]
            a3 = coeffs[3]
            a4 = coeffs[4]
            a6 = coeffs[5]

            b2 = a1^2 + 4*a2
            b4 = 2*a4 + a1*a3
            b6 = a3^2 + 4*a6
            b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2

            d = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
            
            if d != 0
                if d == -16*(4*a4^3 + 27*a6^2)
                    E = new{T}()
                    E.short = true
                    E.base_field = AbstractAlgebra.parent(a4)
                    E.a_invars = (a4, a6)
                    E.b_invars = (b4, b6, b8)
                    E.coeff = [ deepcopy(z) for z ∈ coeffs]
                elseif d == (9 - a6)*(1 + 4*a2)^3 - 432*a6^2
                    E = new{T}()
                    E.simp = true
                    E.base_field = AbstractAlgebra.parent(a1)
                    E.a_invars = (a1, a2, a6)
                    E.b_invars = (b2, b6, b8)
                    E.coeff = [ deepcopy(z) for z ∈ coeffs]
                else
                    E = new{T}()
                    E.base_field = AbstractAlgebra.parent(a1)
                    E.a_invars = (a1, a2, a3, a4, a6)
                    E.b_invars = (b2, b4, b6, b8)
                    E.coeff = [ deepcopy(z) for z ∈ coeffs]
                    E.disc = d
                end
            else
                error("Discriminant is zero")
            end
            return E
        end
    end
end

mutable struct EllCrvPt{T}
    degree::Int
    coord::Tuple{Vararg{T, N} where N}
    isinfinite::Bool
    parent::EllCrv{T}

    function EllCrvPt{T}(E::EllCrv{T}, coords::Tuple{Vararg{T, N} where N}, 
        check::Bool = true) where T
        if check
            if ison_curve(E, coords)
                if length(coords) < 2
                    error("Point must have at least two coordinates.")
                elseif length(coords) == 2
                    P = new{T}(coords[1], coords[2], false, E)
                    return P
                elseif length(coords) == 3
                    P = new{T}(coords[1], coords[2], coords[3], false, E)
                    return P
                end
            else
                error("Point is not on curve.")
            end
        else
            P = new{T}(coords[1]..., false, E)
            return P
        end
    end
    
    # Point at infinity
    function EllCrvPt{T}(E::EllCrv{T}) where T
        z = new{T}()
        z.parent = E
        z.isinfinite = true
        return z
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

function (E::EllCrv{T})(coords::Tuple{Vararg{S, N} where N}, check::Bool = true) where 
    {S, T}
    if length(coords) < 2
        error("Point must have at least two coordinates.")
    elseif length(coords) == 2
        print("Point is affine.")
    elseif length(coords) == 3
        print("Point is projective.")
    end

    if S == T 
        AbstractAlgebra.parent(coords[1]) != base_field(E) &&
            error("The elliptic curve and point must be defined over the same
            field.")
            return EllCrvPt{T}(E, coords, check)
    else
        return EllCrvPt{T}(E, map(base_field(E), coords), check)
    end
end

################################################################################
#
#  Field access/ Quotients
#
################################################################################

@inline function parent(P::EllCrvPt) return P.parent end
@inline function isinfinite(P::EllCrvPt) return P.isinfinite end
@inline function isshort(E::EllCrv) return E.short end
@inline function issimp(E::EllCrv) return E.simp end
@inline function parenttype(::Type{AbstractAlgebra.GFElem{Int64}}) end

@inline function base_field(E::EllCrv{T}) where T
    return E.base_field::parenttype(T)
end

@inline function a_invars(E::EllCrv)
    if isdefined(E, :a_invars)
        return [ deepcopy(z) for z ∈ E.a_invars ]
    else
        t = (E.coeff[1], E.coeff[2], E.coeff[3], E.coeff[4], E.coeff[5])
        E.a_invars = t
        return t
    end
end

@inline function b_invars(E::EllCrv)
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

        t = (b2, b4, b6, b8)
        E.b_invars = t
        return t
    end
end

################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, E::EllCrv)
    if E.short
        print(io, "Elliptic curve with equation y^2 = x^3")
        if !iszero(E.coeff[4])
            print(io, " + $(E.coeff[4])*x")
        end
        if !iszero(E.coeff[5])
            print(io, " + $(E.coeff[5])")
        end
        print(io, "\n")
    if E.simp
        print(io, "Elliptic curve with equation y^2 + x*y = x^3")
        if !iszero(E.coeff[2])
            print(io, " + $(E.coeff[2])*x^2")
        end
        if !iszero(E.coeff[5])
            print(io, " + $(E.coeff[5])")
        end
        print(io, "\n")
    else
        print(io, "Elliptic Curve with equation y^2")
        if !iszero(E.coeff[1])
            print(io, " + $(E.coeff[1])*x*y")
        end
        if !iszero(E.coeff[3])
            print(io, " + $(E.coeff[3])*y")
        end
        print(io, " = x^3")
        if !iszero(E.coeff[2])
            print(io, " + $(E.coeff[2])*x^2")
        end
        if !iszero(E.coeff[4])
            print(io, " + $(E.coeff[4])*x")
        end
        if !iszero(E.coeff[5])
            print(io, " + $(E.coeff[5])")
        end
        print(io, "\n")
    end
end

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

function ison_curve(E::EllCrv{T}, coords::Tuple{Vararg{T, N} where N}) where T 
    length(coords) != 2 && error("Point must have at least two coordinates.")

    (a1, a2, a3, a4, a6) = a_invars(E)
    
    if length(coords) == 2
        x = coords[1]
        y = coords[2]
        
        if y^2 + a1*x*y + a3*y == x^3 + a2*x^2 + a4*x + a6
            return true
        else
            return false
        end
    elseif length(coords) == 3
        X = coords[1]
        Y = coords[2]
        Z = coords[3]

        if Y^2*Z + a1*X*Y*Z + a3*Y*Z^2 == X^3 + a2*X^2*Z + a4*X*Z^2 + a6*Z^3
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

    (b2, b4, b6, b8) = b_invars(E)

    d = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
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

    (b2, b4) = b_invars(E)

    c4 = b2^2 - 24*b4

    j = AbstractAlgebra.divexact(c4^3, disc(E))
    E.j = j
    return j::T
end

################################################################################
#
#  Inversion
#
################################################################################

function ⊟(P::EllCrvPt{T}) where T
    E = parent(P)
    coords = P.coord
    if isinfinite(P)
        O = infinity(E)
        return O
    end

    x = coords[1]
    y = -coords[2]

    t = (x, y)

    Q = E(t, false)

    return Q
end

################################################################################
#
#  Addition of Points (Group Law)
#
################################################################################

# Add two points P, Q with projective coordinates
# P := [P[1]:P[2]:P[3]], Q := [Q[1]:Q[2]:Q[3]] or 
# affine coordinates P:= (P[1], P[2]), Q := (Q[1], Q[2])

function ⊞(P::EllCrvPt{T}, Q::EllCrvPt{T}) where T 
    parent(P) != parent(Q) && error("Points must live on the same curve.")
    coords = P.coord
    E = parent(P)
    
    if length(coords) == 2
        if isshort(E)

            (a4, a6) = a_invars(E)
        
        # Is either P or Q the point at infinity?
            if P.isinfinite
                return Q
            elseif Q.isinfinite
                return P
            end

            if P.coord[1] != Q.coord[1]
                m = AbstractAlgebra.divexact(Q.coord[2] - P.coord[2], 
                Q.coord[1] - P.coord[1])
            
                x = m^2 - P.coord[1] - Q.coord[1]
                y = m*(P.coord[1] - x) - P.coord[2]

                t = (x, y)
            elseif P.coord[2] != Q.coord[2]
                return infinity(E)
            elseif P.coord[2] != 0
                m = AbstractAlgebra.divexact(3*(P.coord[1])^2 + a4, 
                2*(P.coord[2]))
            
                x = m^2 - 2*(P.coord[1])
                y = m*(P.coord[1] - x) - P.coord[2]
            
                t = (x, y)
            else
                return infinity(E)
            end
        Erg = E(t, false)
        else
            error("For short form.")
        end
    end
    return Erg
end

################################################################################
#
#  Window Non-adjacent Form (for Scalar Mult.)
#
################################################################################

function NAF_dw(d::Int, w::Int)
    l = 0
    k = Int[]

    if d > 0 && w >= 2
        while d >= 1
            if mod(d, 2) != 0
                k[l] = 2 - mod(d, 2^w) # TODO: signed modulus mods
                d = d - k[l]
            else
                k[l] = 0
            end
            d = div(d, 2)
            l += 1
        end
        # TODO: concatenate digits of array k into WNAF form for d 
        return k
    end
end

################################################################################
#
#  Scalar multiplication
#
################################################################################

function ×(n::Int, P::EllCrvPt)
    E = parent(P)
    coords = P.coord
    O = infinity(E)
    F = base_field(E)

    if issimp(E) == true
        if F <: AbstractAlgebra.GFField && 
            AbstractAlgebra.characteristic(F) == 2
            
            # Efficient affine point quintupling
            # for simple Weierstrass form Elliptic 
            # Curves over binary fields.
            # https://eprint.iacr.org/2017/840.pdf
            if n == 5
                xP = coords[1]
                yP = coords[2]

                (a1, a2, a6) = a_invars(E)

                A = a2
                B = a6

                α = xP^4 + xP^3 + B
                β = α^2 + xP^2*(xP^4 + B)
                γ = α^2*(xP^4 + B) + xP^3*β

                ρ = AbstractAlgebra.divexact(xP^3*β, γ)
                
                x5P = xP + ρ + ρ^2
                y5P = yP + xP + (x5P + xP)*(ρ + xP^2 + A) + 
                AbstractAlgebra.divexact(xP*β*α^2*(β + 
                (xP^4 + B)*(xP^4 + B + yP^2 + xP^2)), γ^2)

                FiveP = E((x5P, y5P), false)
                
                return FiveP
            end
        end
    end

    if n >= 0
        a = n
    else
        a = -n
    end

    # Double and add
    while a != 0
      if mod(a,2) == 0
        a = div(a,2)
        P = P ⊞ P
      else
        a = a - 1
        O = O ⊞ P
      end
    end
  
    if n < 0
      P = ⊟(P)
    end
  
    return O
end

################################################################################
#
#  Division Polynomials (for SEA)
#
################################################################################

# TODO: Make this dynamic
# https://eprint.iacr.org/2002/109.pdf
function division_polynomialE(E::EllCrv, n::Int, x = nothing, y = nothing)
    
    if isshort(E) == false
        error("Must be in short form.")
    end
    
    (a1, a2, a3, a4, a6) = a_invars(E)

    A = numerator(a4)
    B = numerator(a6)
    
    if x === nothing
        Z = AbstractAlgebra.ZZ
        Zx, _x = AbstractAlgebra.PolynomialRing(Z, "x")
        Zxy, y = AbstractAlgebra.PolynomialRing(Zx, "y")
    else
        Zxy = AbstractAlgebra.parent(x)
    end

    if mod(n, 2) == 1
        m = n 
        return AbstractAlgebra.divexact(division_polynomialE(E, m, x, y), 
        division_polynomialE(E, 2, x, y))
    end
end