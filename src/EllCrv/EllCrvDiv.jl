################################################################################
#
#          EllCrv/EllCrvDiv.jl : Divisors on Elliptic Curves and Logarithms of 
#                                   Divisors
#
################################################################################

################################################################################
#
#  Imports
#
################################################################################

import AbstractAlgebra

using Markdown

include("FieldsRings.jl")
include("EllCrv.jl")
include("EllCrvMap.jl")
include("EllCrvModel.jl")

################################################################################
#
#  Exports
#
################################################################################

export EllCrvDivisor, Logarithm

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

# Generic struct for (Weil) divisors
mutable struct EllCrvDivisor{T, P1,...,Ps}
    coeff::Array{T, 1}
    points::Tuple{P1::EllCrvPt,...,Ps::EllCrvPt}
    func_field_without_zero::FuncField
    rat_func::AbstractAlgebra.FieldElem ∈ func_field_without_zero
    degree::Int
    
    # if assoc(rat_func::AbstractAlgebra.FieldElem, coeff::Array{T, 1}) 
    # --> E::EllCrvDivisor
    is_associated::Bool
    
    # blowup(i::AbstractVarietyMap) --> AbstractVariety
    # if exceptionalDiv(blowup) --> EllCrvDivisor
    is_exceptional::Bool
    
    # if val(coeff::Array{T, 1}) >= 0 
    is_effective::Bool

    function EllCrvDivisor{T, P1,...,Ps}(coeffs::Array{T, 1}, 
        points::Tuple{P1::EllCrvPt,...,Ps::EllCrvPt}, 
        check::Bool = true) where {T, P1,...,Ps}
        if check
            if assoc(rat, coeffs)
                ECD = new{T, P1,...,Ps}()
                ECD.is_associated = true
                return ECD
            elseif val(coeffs) >= 0
                ECD = new{T, P1,...,Ps}()
                ECD.is_effective = true
                return ECD
            else
                ECD = new{T, P1,...,Ps}()
                return ECD
            end
        end
    end
end


@doc Markdown.doc"""
    mutable struct Logarithm

An abstract type extending the notion of logarithms for 
divisors on an elliptic curve from that of field elements. Let 
N = |E(Fq)|. Then... 
"""
# Logarithms of Divisors
mutable struct Logarithm
    base_ec_model::EllCrvModel
    base_divisor::EllCrvDivisor
    #...
end

################################################################################
#
#  Constructors for users
#
################################################################################

function EllCurveDiv{T, P1,...,Ps}(coeff::Array{T}, 
    points::Tuple{P1::EllCrvPt,...,Ps::EllCrvPt}, check::Bool = true) where 
    {T, P1,...,Ps}
        EDiv = EllCrvDivisor{T, P1,...,Ps}(coeff, points, check)
        return EDiv
    end
end

# blowup(i::AbstractVarietyMap) --> V::AbstractVariety

################################################################################
#
#  Functions for Divisor types
#
################################################################################

function ord{P1,...,Ps}(rat_func::AbstractAlgebra.FieldElem, 
    points::Tuple{P1::EllCrvPt,...,Ps::EllCrvPt}) where 
    {P1,...,Ps}
    ##...
    return
end

function assoc{T}(rat_func::AbstractAlgebra.FieldElem, coeff::Array{T, 1}) where T
    coeff = ord(rat_func)
    ##...

    return EllCrvDivisor(coeff, points, true)
end



