################################################################################
#
#          EllCrv/FieldsRings.jl : Fields of definition for Elliptic Curves,
#           Field Extensions --> Number Fields
#           Local Rings --> Function Fields
#
################################################################################


################################################################################
#
#  Submodules
#
################################################################################

include("..src/EllCrv/EllCrv.jl")

################################################################################
#
#  Exports
#
################################################################################

export ExtField, NumField, Loc, LocElem, FuncField

################################################################################
#
#  Abstract Types/Structs
#
################################################################################

abstract type AbstractField{T} <: AbstractAlgebra.Field end
mutable struct ExtField{T} <: AbstractField{T}
    degree::Int
    base_field::AbstractAlgebra.Field
end

mutable struct NumField{T} <: AbstractAlgebra.SimpleNumField{T}
   base_field::AbstractAlgebra.Field
   pol::AbstractAlgebra.Poly{T}
   S::Symbol
   primitive::Bool
   primitive_element::AbstractAlgebra.NumFieldElem

   ###
end

mutable struct Loc{T} <: AbstractAlgebra.Ring
    base_ring::AbstractAlgebra.Ring
    prime::T
    primes::Array{T,1}  # in general, not set.
    comp::Bool  # false: den has to be coprime to prime
                # true:  den can ONLY use prime (and powers)
 
    function Loc{T}(prime::T, primes::Array{T,1}, cached::Bool = true, comp::Bool = false) where 
     {T <: AbstractAlgebra.RingElem}
       length(primes) == 0 && error("No element to localize at since array of primes is empty")
       if cached && haskey(LocDict, (AbstractAlgebra.parent(prime), prime, comp))
          return LocDict[AbstractAlgebra.parent(prime), prime, comp]::Loc{T}
       else
          z = new(AbstractAlgebra.parent(prime), prime, primes, comp)
          if cached
             LocDict[AbstractAlgebra.parent(prime), prime, comp] = z
          end
          return z
       end
    end
 
    function Loc{T}(prime::T, cached::Bool = true, comp::Bool = false) where 
     {T <: AbstractAlgebra.RingElem}
      isunit(prime) && error("no-point")
      if cached && haskey(LocDict, (AbstractAlgebra.parent(prime), prime, comp))
        return LocDict[AbstractAlgebra.parent(prime), prime, comp]::Loc{T}
      else
        r = new()
        r.base_ring = AbstractAlgebra.parent(prime)
        r.prime = prime
        r.comp = comp
        if cached
           LocDict[AbstractAlgebra.parent(prime), prime, comp] = r
        end
        return r
      end
    end
 end

const LocDict = Dict{Tuple{AbstractAlgebra.Ring, AbstractAlgebra.RingElement, Bool}, 
AbstractAlgebra.Ring}()

function isin(a, L::Loc{T}) where {T <: AbstractAlgebra.RingElem}
  iszero(a) && return true
  L.comp || (!isone(gcd(denominator(a), prime(L))) && return false)
  L.comp && ppio(denominator(a), prime(L))[1] != denominator(a.data) && return false
  return true
end

mutable struct LocElem{T} <: AbstractAlgebra.RingElem
   data::AbstractAlgebra.FieldElem
   parent::Loc{T}
   function LocElem{T}(data::AbstractAlgebra.FracElem{T}, par::Loc, checked::Bool = true) where 
    {T <: AbstractAlgebra.RingElem}
     checked && (isin(data, par) || error("illegal elt"))
     return new{T}(data,par)
   end
end

# Function fields of curves.
mutable struct FuncField{T} <: Loc{T}
    # The function field K[E] of the elliptic curve E 
    # consists of ratios of polynomials. Taking an 
    # elliptic curve E as an integral, regular scheme of 
    # finite type, then the function field K[E] is 
    # just the local ring of the generic point.
    
    base_curve::EllCrv
    generic_point::EllCrvPt

    function FuncField{T}()
end

###############################################################################
#
#   Data type and parent object methods
#
###############################################################################

elem_type(::Type{Loc{T}}) where {T} = LocElem{T}

parent_type(::Type{LocElem{T}}) where {T} = Loc{T}

base_ring(L::Loc) = L.base_ring

base_ring(a::LocElem) = base_ring(parent(a))

parent(a::LocElem) = a.parent

function check_parent(a::LocElem{T}, b::LocElem{T})  where 
    {T <: AbstractAlgebra.RingElem}
    parent(a) !== parent(b) && error("Parent objects do not match")
end