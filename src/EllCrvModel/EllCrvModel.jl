import AbstractAlgebra

include("../EllCrv/EllCrv.jl")

export EllCrvModel

mutable struct EllCrvModel{T, P_1,...,P_s}
    elliptic_curve::EllCrv{T}
    ec_point::EllCrvPt{T}
    divisor::EllCrvDivisor{T, P_1,...,P_s}
end
