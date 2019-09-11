export EllCrv, EllCrvPt

export base_field, discriminant, division_polynomial, EllipticCurve, infinity,
    isinifinite, isshort, isoncurve, j_invariant, Psi_polynomial,
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

    torsion_points::Array{EllCrvPt, 1}
    torsion_structure::Tuple{Array{Int, 1}, Array{EllCrvPt, 1}}

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
    coordx::T 
    coordy::T
    isinifinite::Bool
    parent::EllCrv{T}

    function EllCrvPt{T}(E::EllCrv{T}, coords::Array{T, 1}, check::Bool = true) where {T}
        if check
            if ison_curve(E, coords)
                P = new{T}(coords[1], coords[2], false, E)
                return P 
            else
                error("Point is not on curve.")
            end
        else
            P = new{T}(coords[1], coords[2], false, E)
            return P 
        end
    end

    function EllCrvPt{T}(E::EllCrv{T}) where {T}
        z = new{T}()
        z.parent = E 
        z.isinifinite = true
        return z
    end
end