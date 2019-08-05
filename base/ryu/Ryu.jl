module Ryu

include("utils.jl")
include("shortest.jl")
include("fixed.jl")
include("exp.jl")

function writeshortest(x::T) where {T <: Base.IEEEFloat}
    buf, pos = writeshortest(x, zeros(UInt8, 25), 1)
    return String(buf[1:pos-1])
end

function writefixed(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writefixed(x, precision, zeros(UInt8, 2000), 1)
    return String(buf[1:pos-1])
end

function writeexp(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writeexp(x, precision, zeros(UInt8, 2000), 1)
    return String(buf[1:pos-1])
end


end # module