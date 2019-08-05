@inline function writeshortest(x::T, buf::Vector{UInt8}, pos) where {T}
    neg = signbit(x)
    # special cases
    @inbounds if x == 0
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('0')
        buf[pos + neg + 1] = UInt8('e')
        buf[pos + neg + 2] = UInt8('0')
        return buf, pos + neg + 3
    elseif isnan(x)
        buf[pos] = UInt8('N')
        buf[pos + 1] = UInt8('a')
        buf[pos + 2] = UInt8('N')
        return buf, pos + 3
    elseif !isfinite(x)
        if neg
            buf[pos] = UInt8('-')
        end
        buf[pos + neg] = UInt8('I')
        buf[pos + neg + 1] = UInt8('n')
        buf[pos + neg + 2] = UInt8('f')
        return buf, pos + neg + 3
    end

    bits = uint(x)
    mant = bits & (oftype(bits, 1) << mantissabits(T) - oftype(bits, 1))
    exp = Int((bits >> mantissabits(T)) & ((1 << exponentbits(T)) - 1))

    m2 = oftype(bits, 1 << mantissabits(T)) | mant
    e2 = exp - bias(T) - mantissabits(T)
    fraction = m2 & ((oftype(bits, 1) << -e2) - 1)
    if e2 > 0 || e2 < -52 || fraction != 0
        if exp == 0
            e2 = 1 - bias(T) - mantissabits(T) - 2
            m2 = mant
        else
            e2 -= 2
        end
        even = (m2 & 1) == 0
        mv = oftype(m2, 4 * m2)
        mp = oftype(m2, mv + 2)
        mmShift = mant != 0 || exp <= 1
        mm = oftype(m2, mv - 1 - mmShift)
        vmIsTrailingZeros = false
        vrIsTrailingZeros = false
        lastRemovedDigit = 0x00
        if e2 >= 0
            q = log10pow2(e2) - (T == Float64 ? (e2 > 3) : 0)
            e10 = q
            k = pow5_inv_bitcount(T) + pow5bits(q) - 1
            i = -e2 + q + k
            vr, vp, vm = mulshiftinvsplit(T, mv, mp, mm, q, i)
            if T == Float32 || T == Float16
                if q != 0 && div(vp - 1, 10) <= div(vm, 10)
                    l = pow5_inv_bitcount(T) + pow5bits(q - 1) - 1
                    @inbounds mul = T == Float32 ? FLOAT_POW5_INV_SPLIT[q] : HALF_POW5_INV_SPLIT[q]
                    lastRemovedDigit = (mulshift(mv, mul, -e2 + q - 1 + l) % 10) % UInt8
                end
            end
            if q <= qinvbound(T)
                if ((mv % UInt32) - 5 * div(mv, 5)) == 0
                    vrIsTrailingZeros = pow5(mv, q)
                elseif even
                    vmIsTrailingZeros = pow5(mm, q)
                else
                    vp -= pow5(mp, q)
                end
            end
        else
            q = log10pow5(-e2) - (T == Float64 ? (-e2 > 1) : 0)
            e10 = q + e2
            i = -e2 - q
            k = pow5bits(i) - pow5_bitcount(T)
            j = q - k
            vr, vp, vm = mulshiftsplit(T, mv, mp, mm, i, j)
            if T == Float32 || T == Float16
                if q != 0 && div(vp - 1, 10) <= div(vm, 10)
                    j = q - 1 - (pow5bits(i + 1) - pow5_bitcount(T))
                    @inbounds mul = T == Float32 ? FLOAT_POW5_SPLIT[i + 2] : HALF_POW5_SPLIT[i + 2]
                    lastRemovedDigit = (mulshift(mv, mul, j) % 10) % UInt8
                end
            end
            if q <= 1
                vrIsTrailingZeros = true
                if even
                    vmIsTrailingZeros = mmShift
                else
                    vp -= 1
                end
            elseif q < qbound(T)
                vrIsTrailingZeros = pow2(mv, q - (T != Float64))
            end
        end
        removed = 0
        if vmIsTrailingZeros || vrIsTrailingZeros
            while true
                vpDiv10 = div(vp, 10)
                vmDiv10 = div(vm, 10)
                vpDiv10 <= vmDiv10 && break
                vmMod10 = (vm % UInt32) - UInt32(10) * (vmDiv10 % UInt32)
                vrDiv10 = div(vr, 10)
                vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                vmIsTrailingZeros &= vmMod10 == 0
                vrIsTrailingZeros &= lastRemovedDigit == 0
                lastRemovedDigit = vrMod10 % UInt8
                vr = vrDiv10
                vp = vpDiv10
                vm = vmDiv10
                removed += 1
            end
            if vmIsTrailingZeros
                while true
                    vmDiv10 = div(vm, 10)
                    vmMod10 = (vm % UInt32) - UInt32(10) * (vmDiv10 % UInt32)
                    vmMod10 != 0 && break
                    vpDiv10 = div(vp, 10)
                    vrDiv10 = div(vr, 10)
                    vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                    vrIsTrailingZeros &= lastRemovedDigit == 0
                    lastRemovedDigit = vrMod10 % UInt8
                    vr = vrDiv10
                    vp = vpDiv10
                    vm = vmDiv10
                    removed += 1
                end
            end
            if vrIsTrailingZeros && lastRemovedDigit == 5 && vr % 2 == 0
                lastRemovedDigit = UInt8(4)
            end
            output = vr + ((vr == vm && (!even || !vmIsTrailingZeros)) || lastRemovedDigit >= 5)
        else
            roundUp = false
            vpDiv100 = div(vp, 100)
            vmDiv100 = div(vm, 100)
            if vpDiv100 > vmDiv100
                vrDiv100 = div(vr, 100)
                vrMod100 = (vr % UInt32) - UInt32(100) * (vrDiv100 % UInt32)
                roundUp = vrMod100 >= 50
                vr = vrDiv100
                vp = vpDiv100
                vm = vmDiv100
                removed += 2
            end
            while true
                vpDiv10 = div(vp, 10)
                vmDiv10 = div(vm, 10)
                vpDiv10 <= vmDiv10 && break
                vrDiv10 = div(vr, 10)
                vrMod10 = (vr % UInt32) - UInt32(10) * (vrDiv10 % UInt32)
                roundUp = vrMod10 >= 5
                vr = vrDiv10
                vp = vpDiv10
                vm = vmDiv10
                removed += 1
            end
            output = vr + (vr == vm || roundUp || lastRemovedDigit >= 5)
        end
        nexp = e10 + removed
    else
        output = m2 >> -e2
        nexp = 0
        while true
            q = div(output, 10)
            r = (output % UInt32) - UInt32(10) * (q % UInt32)
            r != 0 && break
            output = q
            nexp += 1
        end
    end

    if neg
        @inbounds buf[pos] = UInt8('-')
        pos += 1
    end
    olength = decimallength(output)

    i = 0
    if (output >> 32) != 0
        q = output ÷ 100000000
        output2 = (output % UInt32) - UInt32(100000000) * (q % UInt32)
        output = q

        c = output2 % 10000
        output2 = div(output2, 10000)
        d = output2 % 10000
        c0 = (c % 100) << 1
        c1 = (c ÷ 100) << 1
        d0 = (d % 100) << 1
        d1 = (d ÷ 100) << 1
        unsafe_copyto!(buf, pos + olength - 1, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + olength - 3, DIGIT_TABLE, c1 + 1, 2)
        unsafe_copyto!(buf, pos + olength - 5, DIGIT_TABLE, d0 + 1, 2)
        unsafe_copyto!(buf, pos + olength - 7, DIGIT_TABLE, d1 + 1, 2)
        i += 8
    end
    output2 = output % UInt32
    while output2 >= 10000
        c = output2 % 10000
        output2 = div(output2, 10000)
        c0 = (c % 100) << 1
        c1 = (c ÷ 100) << 1
        unsafe_copyto!(buf, pos + olength - i - 1, DIGIT_TABLE, c0 + 1, 2)
        unsafe_copyto!(buf, pos + olength - i - 3, DIGIT_TABLE, c1 + 1, 2)
        i += 4
    end
    if output2 >= 100
        c = (output2 % 100) << 1
        output2 = div(output2, 100)
        unsafe_copyto!(buf, pos + olength - i - 1, DIGIT_TABLE, c + 1, 2)
        i += 2
    end
    if output2 >= 10
        c = output2 << 1
        #=@inbounds=# buf[pos + olength - i] = DIGIT_TABLE[c + 2]
        #=@inbounds=# buf[pos] = DIGIT_TABLE[c + 1]
    else
        #=@inbounds=# buf[pos] = UInt8('0') + (output2 % UInt8)
    end

    if olength > 1
        #=@inbounds=# buf[pos + 1] = UInt8('.')
        pos += olength + 1
    else
        pos += 1
    end

    #=@inbounds=# buf[pos] = UInt8('e')
    pos += 1
    exp2 = nexp + olength - 1
    if exp2 < 0
        #=@inbounds=# buf[pos] = UInt8('-')
        pos += 1
        exp2 = -exp2
    end

    if exp2 >= 100
        c = exp2 % 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * div(exp2, 10) + 1, 2)
        #=@inbounds=# buf[pos + 2] = UInt8('0') + (c % UInt8)
        pos += 3
    elseif exp2 >= 10
        unsafe_copyto!(buf, pos, DIGIT_TABLE, 2 * exp2 + 1, 2)
        pos += 2
    else
        #=@inbounds=# buf[pos] = UInt8('0') + (exp2 % UInt8)
        pos += 1
    end

    return buf, pos
end