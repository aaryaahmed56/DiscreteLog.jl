import AbstractAlgebra

export EllCrv, EllCrvPt, EllCrvDivisor

export base_field, discriminant, division_polynomial, EllipticCurve, infinity,
    is_inifinite, isshort, is_on_curve, j_invariant, Psi_polynomial, is_rational,
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
    is_inifinite::Bool
    parent::EllCrv{T}

    function EllCrvPt{T}(E::EllCrv{T}, coords::Array{T, 1}, check::Bool = true) where {T}
        if check
            if is_on_curve(E, coords)
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

    function EllCrvDivisor{T, P_1,...,P_s}(rat::AbstractAlgebra.FieldElem, coeffs::Array{T, 1}, 
        check::Bool = true) where {T, P_1,...,P_s}
        if check
            if is_rational(rat, coeffs)
                ECD = new(){T, P_1,...,P_s}
                ECD.is_associated = true
            end
        end
    end
end
