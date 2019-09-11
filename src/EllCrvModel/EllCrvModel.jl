import AbstractAlgebra

include("../EllCrv/EllCrv.jl")

export EllCrvModel

mutable struct EllCrvModel{T, T, T}
    elliptic_curve::EllCrv{T}
    ec_point::EllCrvPt{T}
    divisor::T
end
