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

include("EllCrv.jl")
include("EllCrvModel.jl")

################################################################################
#
#  Exports
#
################################################################################

export EllCrvDivisor, Logarithm

mutable struct EllCrvDivisor{T, P_1,...,P_s}
    coeff::Array{T, 1}
    s = length(coeff)
    points::Tuple{P_1::EllCrvPt{T},...,P_s::EllCrvPt{T}}
    func_field_without_zero::AbstractAlgebra.Field
    rat_func::AbstractAlgebra.FieldElem ∈ func_field_without_zero
    degree::Int
    is_associated::Bool
    
    # Effective divisors are imporant to track as the elimination
    # procedures presume we have a divisor D on the elliptic curve 
    # that is not part of a collection of effective divisors.
    is_effective::Bool
    
    is_exceptional::Bool

    function EllCrvDivisor{T, P_1,...,P_s}(rat::AbstractAlgebra.FieldElem, 
        coeffs::Array{T, 1}, points::Tuple{P_1::EllCrvPt{T},...,P_s::EllCrvPt{T}}, 
        check::Bool = true) where {T, P_1,...,P_s}
        if check
            # if is_rational(rat, coeffs)
                ECD = new{T, P_1,...,P_s}()
                ECD.is_associated = true
                return ECD
        else 
            ECD = new{T, P_1,...,P_s}()
            ECD.is_exceptional = true
            return ECD
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

#function EllipticCurveDiv(x::Array{T}, y::Tuple{P_1,...,P_s}, check::Bool = true) 
#    where {P_1,...,P_s, T}
#    EDiv = EllCrvDivisor{T, P_1,...,P_s, check}
#    return EDiv
#end