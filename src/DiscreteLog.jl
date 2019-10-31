__precompile__()

module DiscreteLog

################################################################################
#
#  Import/export
#
################################################################################

using LinearAlgebra, Distributed, SparseArrays, Serialization, Random, Pkg, Test

import AbstractAlgebra

import SparseArrays

import Serialization: serialize, deserialize

###############################################################################
#
#   Library initialization
#
###############################################################################

const pkgdir = joinpath(dirname(pathof(DiscreteLog)), "..")

function __init__()

    if myid() == 1
        println("")
        print("Welcome to DiscreteLog.jl")
        println()
    end
end

################################################################################
#
#  Submodules
#
################################################################################

include("../src/ElimProcedures/Deg32Elim.jl")
include("../src/ElimProcedures/Deg43Elim.jl")

include("../src/EllCrv/EllCrv.jl")
include("../src/EllCrv/EllCrvDiv.jl")
include("../src/EllCrv/EllCrvFin.jl")
include("../src/EllCrv/EllCrvMap.jl")

include("../src/Misc/BigIntInterface.jl")
include("../src/Misc/FieldsRings.jl")
include("../src/Misc/PolyFac.jl")

end