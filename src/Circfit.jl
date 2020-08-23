module Circfit
using LsqFit

export circfit

"""
x: Array of x-coordinates
y: Array of y-coordinates
p0: Array of initial values [center_x, center_y, radius] 
"""
function circfit(x::AbstractArray,y::AbstractArray,p0::AbstractArray,kwargs...)
    T = eltype(y)
    f = (p) -> @. p[3]^2 - (x-p[1])^2 - (y-p[2])^2
    LsqFit.lmfit(f,p0,T[];kwargs...)
end

end # module
