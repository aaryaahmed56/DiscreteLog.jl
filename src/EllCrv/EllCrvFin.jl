################################################################################
#
#          EllCrv/EllCrvFin.jl : Elliptic curves over finite fields
#
#
# Copyright (c) 2015, 2016: Claus Fieker, Tommy Hofmann
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# (C) 2016 Tommy Hofmann
# (C) 2016 Robin Ammon
# (C) 2016 Sofia Brenner
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
using Random

################################################################################
#
#  Exports
#
################################################################################

export order, order_via_schoof, rand_point, crypto_curve

################################################################################
#
#  Random Point Generator
#
################################################################################

function rand_point(E::EllCrv)
    k = base_field(E)

    if E.short == false
        error("Not for long form.")
    end

    while true

        x = rand(k)
        square = x^3 + E.coeff[4]*x + E.coeff[5]

        a = AbstractAlgebra.issquare(square)
        if a[1] == true
            y = a[2]
            P = E((x, y))
            return P
        end
    end
end

################################################################################
#
#  Elliptic Curve Generator over Finite Fields
#  https://pdfs.semanticscholar.org/c60c/46a9ee6f90958f609f464723f6bf92835feb.pdf
#
################################################################################

function crypto_curve(r0::Union{Int128, BigInt}, k0::Int, h0::Int)
    if prod(r0, k0, h0) <= 0
        error("r0, k0, h0 must be positive integers.")
    end

    if r0 < 2^BigInt(159) || k0 < 4 || h0 < 200
        error("Invalid input.")
    end

    (Δ, p, r, k) = find_prime(r0, k0, h0)
    
    (E, G) = find_curve(Δ, p, r, k)

    return (Δ, p, r, k, E, G)
end

function find_prime(r0::Union{Int128, BigInt}, k0::Int, h0::Int)
    print("Enter a rational prime p: ")
    input = readline()
    p = parse(T, chomp(input)) where T #TODO: Implement Rational Prime type
    
    m = 50
    
    if p == 0
        find_prime_delta_fixed(r0, k0, h0)
    elseif is_prime(p, m) == true && floor(log2(p)) == floor(log2(r0*k0))
        find_discriminant(p, r0, k0, h0)
    else
        error("Invalid prime.")
    end
end

function find_prime_delta_fixed(r0::Union{Int128, BigInt}, k0::Int, h0::Int)
    if k0 > 4
        b = floor()