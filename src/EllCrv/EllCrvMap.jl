################################################################################
#
#          EllCrv/EllCrvMap.jl : Abstract Variety Maps,
#               Elliptic Curve Maps, and
#               Elliptic Curve Blowups at subvarieties.
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

export AbstractVarietyMap, AbstractCurveMap, EllCrvBlowup, EllCrvMap

################################################################################
#
#  Abstract Types/Structs
#
################################################################################
Map(::Type{T}) where T <: AbstractAlgebra.Map = AbstractAlgebra.Map(T)

mutable struct AbstractVarietyMap <: AbstractAlgebra.Map{AbstractVariety, 
    AbstractVariety}
    source::AbstractVariety
    target::AbstractVariety

    function AbstractVarietyMap(X::AbstractVariety, Y::AbstractVariety)
        z = new(source, target)
        return z
    end
end

mutable struct AbstractCurveMap <: AbstractVarietyMap
    source::AbstractCurve
    target::AbstractCurve

    function AbstractCurveMap(X::AbstractCurve, Y::AbstractCurve)
        z = new(source, target)
        return z
    end
end

mutable struct EllCrvBlowup <: AbstractVarietyMap
    source::AbstractVariety
    target::EllCrv

    function EllCrvBlowup(Z::AbstractVariety, Y::EllCrv)
        z = new(source, target)
        return z
    end
end

mutable struct EllCrvMap <: AbstractCurveMap
    source::EllCrv
    target::EllCrv

    function EllCrvMap(E1::EllCrv, E2::EllCrv)
        z = new(source, target)
        return z
    end
end
