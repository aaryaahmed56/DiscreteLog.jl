import AbstractAlgebra

using Markdown

include("EllCrv.jl")

export EllCrvMap, EllCrvBlowup

mutable struct AbstractVarietyMap{F}
end

mutable struct EllCrvMap{F} <: AbstractVarietyMap
end

mutable struct EllCrvBlowup{F} <: EllCrvMap
end
