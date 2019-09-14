  
################################################################################
#
#          EllCrv/EllCrvDiv.jl : Divisors on Elliptic Curves and Logarithms of 
#                                   Divisors
#
################################################################################


################################################################################
#
#  Submodules
#
################################################################################

include("FieldsRings.jl")
include("EllCrv.jl")

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

export EllCrvDivisor, Logarithm
export EllipticCurveDivisor, ord, assoc

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

# Generic struct for Weil divisors
mutable struct EllCrvDivisor{T, P1,...}
    degree::Int
    coeff::Array{T, 1}
    points::Tuple{P1::EllCrvPt,...}
    func_field_without_zero::FuncField
    rat_func::AbstractAlgebra.FieldElem ∈ func_field_without_zero
    
    # if assoc(rat_func::AbstractAlgebra.FieldElem, coeff::Array{T, 1}) 
    # --> EllCrvDivisor
    is_associated::Bool

    # blowup(i::AbstractVarietyMap) --> AbstractVariety
    # if exceptionalDiv(blowup) --> EllCrvDivisor
    is_exceptional::Bool

    # if cumprod(coeff::Array{T, 1}) >= 0 
    is_effective::Bool

    function EllCrvDivisor{T, P1,...}(coeffs::Array{T, 1}, 
        points::Tuple{P1::EllCrvPt,...}, check::Bool = true) where {T, P1,...}
        if check
            if assoc(rat_func, coeffs, points)
                ECD = new{T, P1,...}()
                ECD = coeffs[1]*P1 +...
                ECD.is_associated = true
            elseif cumprod(coeff) >= 0
                ECD = new{T, P1,...}()
                ECD = coeffs[1]*P1 +...
                ECD.is_effective = true
            end
        else 
            ECD = new{T, P1,...}()
            ECD = coeffs[1]*P1 +...
        end
        return ECD
    end
end

# Logarithms of Divisors
mutable struct Logarithm
    base_ec_divisor::EllCrvDivisor
    #...
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllipticCurveDivisor{T, P1,...}(coeff::Array{T, 1}, 
    points::Tuple{P1::EllCrvPt,...}, check::Bool = true) where {T, P1,...}
    if check
        ECD = EllCrvDivisor{T, P1,...}(coeff, points, check)
        return ECD
    end
end

################################################################################
#
#  Divisorial Types
#
################################################################################

function ord{P1,..}(rat_func::AbstractAlgebra.FieldElem, 
    points::Tuple{P1::EllCrvPt,...}) where {P1,...}
    #...
    return
end

function assoc{T}(rat_func::AbstractAlgebra.FieldElem, 
    coeff::Array{T, 1}, points::Tuple{P1::EllCrvPt,...}) where {T, P1,...}

    coeff = ord(rat_func)
    return EllipticCurveDivisor(coeff, points, true)
end