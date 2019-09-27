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

include("BigIntInterface.jl")
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
    (a4, a6) = a_invars(E)

    if E.short == false
        error("Not implemented for long form.")
    end

    while true

        x = rand(k)
        square = x^3 + a4*x + a6

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
#  Elliptic Curve Generator over Finite Fields (Naive implementation at first.)
#  https://pdfs.semanticscholar.org/c60c/46a9ee6f90958f609f464723f6bf92835feb.pdf
#
################################################################################

function crypto_curve(r0::BigInt, k0::BigInt, h0::BigInt)
    if prod(r0, k0, h0) <= 0
        error("r0, k0, h0 must be positive integers.")
    end

    if r0 < 2^BigInt(159) || k0 < 4 || h0 < 200
        error("Invalid input.")
    end

    (d, p, r, k) = find_prime(r0, k0, h0)
    
    (E, G) = find_curve(d, p, r, k)

    return (d, p, r, k, E, G)
end

function find_prime(r0::BigInt, k0::BigInt, h0::BigInt)
    print("Enter a (rational) prime p: ")
    input = readline()
    p = parse(BigInt, chomp(input))
    
    m = 50
    
    if p == 0
        find_prime_delta_fixed(r0, k0, h0)
    elseif is_prime(p, m) == true && floor(log2(p)) == floor(log2(r0*k0))
        find_discriminant(p, r0, k0, h0)
    else
        error("Invalid prime.")
    end
end

function find_prime_delta_fixed(r0::BigInt, k0::BigInt, h0::BigInt)
    print("Enter a discriminant: ")
    input = readline()
    d = parse(BigInt, chomp(input))

    abs_d = abs(d)
    
    if k0 >= 4
        if mod(d, 8) == 1 || mod(abs_d, 8) == 7
            b = BigInt(floor(log2(4*r0)))
            i0 = 2^b
            i1 = 2^(b + 1) - 1
            interval = collect(BigInt, i0:i1)
            
            if find_prime_1_mod_8(r0, 4, d, interval)
                return true
            else
                error("Unable to find prime.")
            end
        else
            error("Get discriminant congruent to 1 mod 8, 
            or its norm congruent to 7 mod 8.")
        end
    elseif k0 >= 2
        if mod(d, 16) == 8 || 12
            b = BigInt(floor(log2(r0*k0)))
            i0 = 2^b
            i1 = 2^(b + 1) - 1
            interval = collect(BigInt, i0:i1)

            if find_prime_0_mod_4(r0, k0, d, interval)
                return true
            else
                error("Unable to find prime.")
            end
        else
            error("Get discriminant congruent to 8 or 12 modulo 16.")
        end
    else
        if mod(d, 8) == 5
            b = BigInt(floor(log2(r0)))
            i0 = 2^b
            i1 = 2^(b + 1) - 1
            interval = collect(BigInt, i0:i1)

            if find_prime_5_mod_8(r0, 1, d, interval)
                return true
            else
                error("Unable to find prime.")
            end
        else
            error("Get discrimint congruent to 5 mod 8.")
        end
    end
end

function find_prime_0_mod_4(r0::BigInt, k0::BigInt, d::BigInt, 
    interval::Array{BigInt, 1})
    # Number of successive tested pairs (t, y).
    # If not successful after T tries, init a new pair.
    const T = 2000

    # Bit of 2 in binary representation of t.
    t_bit_1 = Int

    abs_delta = SpBigInt(d)
    
    ###
end



