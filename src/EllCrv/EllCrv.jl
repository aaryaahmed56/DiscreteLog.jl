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
export base_field, disc, EllipticCurve, infinity, is_infinite, 
is_on_curve, j_invariant, ⊟, ⊞, *

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

abstract type AbstractVariety{T} end

abstract type AbstractCurve{T} <: AbstractVariety{T} end
mutable struct EllCrv{T} <: AbstractCurve{T}
    base_field::AbstractAlgebra.Field
    coeff::Array{T, 1}
    a_invars::Tuple{Vararg{T, N} where N}
    b_invars::Tuple{Vararg{T, N} where N}
    disc::T
    j::T

    function EllCrv{T}(coeffs::Array{T, 1}, check::Bool = true) where T

        if check
            d = 4*coeffs[1]^3 + 27*coeffs[2]^2
            if d != 0
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
    degree::Int
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
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllipticCurve(coeffs::Array{T, 1}, check::Bool = true) where T
    E = EllCrv{T}(coeffs, check)
    return E
end

function EllipticCurvePoint(E::EllCrv{T}, coords::Array{S, 1}, check::Bool = true) where {S, T}
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
#  Field access/ Quotients
#
################################################################################

@inline function parent_type(::Type{AbstractAlgebra.FieldElem}) end
@inline function parent(P::EllCrvPt) return P.parent end
@inline function is_infinite(P::EllCrvPt) return P.is_infinite end

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
        # Coefficients for general Weierstrass equation.
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
        
        if y^2 == x^3 + (E.coeff[1])*x + (E.coeff[2]) || 
            y^2 + x*y == x^3 + (E.coeff[1])*x^2 + (E.coeff[2])
            return true
        else
            return false
        end
    elseif length(coords) == 3
        X = coords[1]
        Y = coords[2]
        Z = coords[3]

        if Y^2*Z == X^3 + (E.coeff[1])*X*Z^2 + (E.coeff[2])*Z^3
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
#  Inverse
#
################################################################################

function ⊟(P::EllCrvPt{T}) where T
    E = P.parent
    coords = P.coord
    if is_infinite(P)
        O = infinity(E)
        return O
    end

    Q = E([coords[1], -coords[2]], false)

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
    
    if length(coords) == 2
        
        # Is either P or Q the point at infinity?
        if P.is_infinite
            return Q
        elseif Q.is_infinite
            return P
        end
        
        E = P.parent

        if P.coord[1] != Q.coord[1]
            m = AbstractAlgebra.divexact(Q.coord[2] - P.coord[2], 
            Q.coord[1] - P.coord[1])
            
            x = m^2 - P.coord[1] - Q.coord[1]
            y = m*(P.coord[1] - x) - P.coord[2]
        elseif P.coord[2] != Q.coord[2]
            return infinity(E)
        elseif P.coord[2] != 0
            m = AbstractAlgebra.divexact(3*(P.coord[1])^2 + (E.coeff[1]), 
            2*(P.coord[2]))
            
            x = m^2 - 2*(P.coord[1])
            y = m*(P.coord[1] - x) - P.coord[2]
        else
            return infinity(E)
        end
        
        Erg = E([x, y], false)
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
                kl = 2 - mod(d, 2^w)
                push!(k, kl)
                d = d - kl
            else
                kl = 0
                push!(k, kl)
            end
            d = div(d, 2)
            l += 1
        end 
        return k 
    end
end

################################################################################
#
#  Scalar multiplication
#
################################################################################

function *(n::Int, P::EllCrvPt)
    E = P.parent
    coords = P.coord
    O = infinity(E)
    F = base_field(E)

    if length(coords) == 2
        if F <: AbstractAlgebra.GFField && 
            AbstractAlgebra.characteristic(F) == 2
            
            # Efficient affine point quintupling
            # https://eprint.iacr.org/2017/840.pdf
            if n == 5 && 6*P != O
                xP = coords[1]
                yP = coords[2]

                A = E.coeff[1]
                B = E.coeff[2]

                α = xP^4 + xP^3 + B
                β = α^2 + xP^2*(xP^4 + B)
                γ = α^2*(xP^4 + B) + xP^3*β

                ρ = AbstractAlgebra.divexact(xP^3*β, γ)
                
                x5P = xP + ρ + ρ^2
                y5P = yP + xP + (x5P + xP)*(ρ + xP^2 + A) + 
                AbstractAlgebra.divexact(xP*β*α^2*(β + 
                (xP^4 + B)*(xP^4 + B + yP^2 + xP^2)), γ^2)

                FiveP = E([x5P, y5P], false)
                
                return FiveP
            end
        end
    end

    # Repeated addition
    if n >= 0
        a = n
    else
        a = -n
    end

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