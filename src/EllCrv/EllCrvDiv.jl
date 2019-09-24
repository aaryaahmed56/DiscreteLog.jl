  
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

export AbstractDivisor, EllCrvDivisor, EllCrvModel, Logarithm
export EllipticCurveDivisor, ord, assoc

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

abstract type AbstractDivisor{EllCrvPt} end

mutable struct EllCrvDivisor{EllCrvPt} <: AbstractDivisor{EllCrvPt}
    degree::Int
    coeff::Array{Int, 1}
    points::Tuple{Vararg{EllCrvPt, N} where N}
    func_field_without_zero::AbstractAlgebra.Field
    rat_func::AbstractAlgebra.FieldElem
    principal::Bool

    # blowup(i::AbstractVarietyMap) --> AbstractVariety
    # if exceptionalDiv(blowup) --> EllCrvDivisor
    exceptional::Bool

    # if prod(coeff::Array{T, 1}) >= 0 
    effective::Bool

    function EllCrvDivisor{EllCrvPt}(coeffs::Array{Int, 1}, 
        points::Tuple{Vararg{EllCrvPt, N} where N}, check::Bool = true) where 
        EllCrvPt
        if check
            if assoc(rat_func, coeffs, points)
                ECD = new{EllCrvPt}()
                ECD = ⊞(coeffs[1]×EllCrvPt...)
                ECD.principal = true
            elseif prod(coeffs) >= 0
                ECD = new{EllCrvPt}()
                ECD = ⊞(coeffs[1]×EllCrvPt...)
                ECD.effective = true
            end
        else 
            ECD = new{EllCrvPt}()
            ECD = ⊞(coeffs[1]×EllCrvPt...)
        end
        return ECD
    end
end

mutable struct EllCrvModel
    curve::EllCrv
    point::EllCrvPt
    divisor::EllCrvDivisor
end

# Logarithms of Divisors
mutable struct Logarithm
    base_ec_model::EllCrvModel
    divisor::EllCrvDivisor
    #...
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllipticCurveDivisor(coeff::Array{Int, 1}, 
    points::Tuple{Vararg{EllCrvPt, N} where N}, check::Bool = true)
    if check
        ECD = EllCrvDivisor{EllCrvPt}(coeff, points, check)
        return ECD
    end
end

################################################################################
#
#  Divisorial Type Methods
#
################################################################################

function ord(rat_func::AbstractAlgebra.FieldElem, 
    points::Tuple{Vararg{EllCrvPt, N} where N})
    
    if AbstractAlgebra.denominator(rat_func) == 0
        
    end

    return
end

function assoc(rat_func::AbstractAlgebra.FieldElem, 
    coeff::Array{Int, 1}, points::Tuple{Vararg{EllCrvPt, N} where N})

    coeff = ord(rat_func, points)
    return EllipticCurveDivisor(coeff, points, true)
end