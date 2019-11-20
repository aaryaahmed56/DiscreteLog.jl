################################################################################
#
#          Misc/BigIntInterface.jl : BigInt Handling and Conversions
#
################################################################################j


################################################################################
#
#  Inline Bitwise Assignment Methods for BigInt
#
################################################################################

@inline function assign_zero() end::BigInt
@inline function assign_one() end
@inline function assign(Union{::Int64, ::Int128, 
    ::BigInt, ::Float32, ::Float64}) end

################################################################################
#
#  Conversions
#
################################################################################

# Convert BigInt to Hex Array
function int_to_bytes(x::BigInt)
    n_bytes_with_zeros = x.size*sizeof(Sys.WORD_SIZE)
    uint8_ptr = convert(Ptr{UInt8}, x.d)
    n_bytes_without_zeros = 1

    if ENDIAN_BOM == 0x04030201
        for i in n_bytes_with_zeros:-1:1
            if unsafe_load(uint8_ptr, i) != 0x00
                n_bytes_without_zeros = i
                break
            end
        end

        result = Array{UInt8}(undef, n_bytes_without_zeros)

        for i in 1:n_bytes_without_zeros
            @inbounds result[n_bytes_without_zeros + 1 - i] = unsafe_load(uint8_ptr, i)
        end
    else
        for i in 1:n_bytes_with_zeros
            if unsafe_load(uint8_ptr, i) != 0x00
                n_bytes_without_zeros = i
                break
            end
        end

        result = Array{UInt8}(undef, n_bytes_without_zeros)

        for i in 1:n_bytes_without_zeros
            @inbounds result[i] = unsafe_load(uint8_ptr, i)
        end
    end
    return result
end

# Convert Hex Array to BigInt
function bytes_to_big(x::Array{UInt8, 1})
    hex = bytes2hex(x)
    return parse(BigInt, hex, base=16)
end