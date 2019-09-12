__precompile__()

module DiscreteLog

################################################################################
#
#  Import/export
#
################################################################################

using LinearAlgebra, Markdown, InteractiveUtils, Libdl, Distributed, Printf, SparseArrays, Serialization, Random, Pkg, Test

import AbstractAlgebra

import SparseArrays

import Serialization: serialize, deserialize

###############################################################################
#
#   Library initialisation
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
#  "Submodules"
#
################################################################################

include("ElimProcedures/Deg32Elim.jl")
include("ElimProcedures/Deg43Elim.jl")
include("EllCrv/EllCrv.jl")
include("EllCrv/EllCrvDiv.jl")
include("EllCrv/EllCrvModel.jl")
include("EllCrv/FieldsRings.jl")