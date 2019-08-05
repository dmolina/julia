module Ryu

include("utils.jl")
include("shortest.jl")
include("fixed.jl")
include("exp.jl")

shortestdigits(::Type{Float64}) = 25
shortestdigits(::Type{Float32}) = 16
shortestdigits(::Type{Float16}) = 7

function writeshortest(x::T) where {T <: Base.IEEEFloat}
    buf, pos = writeshortest(x, zeros(UInt8, shortestdigits(T)), 1)
    return String(buf[1:pos-1])
end

function writefixed(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writefixed(x, precision, zeros(UInt8, precision + shortestdigits(T)), 1)
    return String(buf[1:pos-1])
end

function writeexp(x::T, precision) where {T <: Base.IEEEFloat}
    buf, pos = writeexp(x, precision, zeros(UInt8, precision + shortestdigits(T)), 1)
    return String(buf[1:pos-1])
end

function Base.show(io::IO, x::T) where {T <: Base.IEEEFloat}
    if get(io, :compact, false)
        precision = T == Float16 ? 5 : 6
        buf, pos = writefixed(x, precision, Vector{UInt8}(undef, precision + shortestdigits(T)), 1)
        GC.@preserve buf unsafe_write(io, pointer(buf), pos - 1)
    else
        buf, pos = writeshortest(x, zeros(UInt8, shortestdigits(T)), 1)
        GC.@preserve buf unsafe_write(io, pointer(buf), pos - 1)
    end
    return
end

end # module