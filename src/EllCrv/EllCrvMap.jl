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

export AbstractVarietyMap, EllCrvMap, EllCrvBlowup

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

mutable struct AbstractVarietyMap
    source, target::Tuple{AbstractVariety, AbstractVariety}
end

mutable struct EllCrvMap <: AbstractVarietyMap
    source, target::Tuple{EllCrv, EllCrv}
end

mutable struct EllCrvBlowup{i, e} <: EllCrvMap
    source, target::Tuple{i::AbstractVariety, e::EllCrv}
end