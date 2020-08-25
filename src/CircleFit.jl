module CircleFit
import Statistics: var, cov, stdm

export circfit

"""
Fit a circle to points provided as arrays of x and y coordinates

Example
```
x = [-1.0,0,0,1]
y = [0.0,1,-1,0]
x0,y0,radius = circfit(x,y)
```
"""
circfit(x,y) = kasa(x,y)

"""
Fit a circle to the points provided as arrays of x and y coordinates

This method uses [Kåsa's method](https://doi.org/10.1109/TIM.1976.6312298)
The result is a GeometryBasics::Circle
"""
function kasa(x::AbstractArray, y::AbstractArray)
    x² = x.^2
    y² = y.^2
    
    A = var(x) 
    B = cov(x, y) 
    C = var(y) 
    D = cov(x, y²) + cov(x, x²)
    E = cov(y, x²) + cov(y, y²) 

    ACB2 = 2 * (A * C - B^2)
    am = (D * C - B * E) / ACB2 
    bm = (A * E - B * D) / ACB2
    rk = hypot(stdm(x, am, corrected=false), stdm(y, bm, corrected=false))

    (am, bm, rk)
end

using LinearAlgebra

"""
Fit a circle by using Taubin's method
https://doi.org/10.1007/s10851-005-0482-8
Warning: not optimized
"""
function taubin(x,y)

    z = x.^2 .+ y.^2
    Mx = sum(x)
    My = sum(y)
    Mz = sum(z)
    Mxx = sum(x.^2)
    Myx = Mxy = sum(x.*y)
    Mzx = Mxz = sum(x.*z)
    Myy = sum(y.^2)
    Mzy = Myz = sum(y.*z)
    Mzz = sum(z.^2)
    n = length(x)

    C = [4Mz 2Mx 2My 0
         2Mx n   0   0
         2My 0   n   0
         0   0   0   0]
        
    M = [Mzz Mxz Myz Mz
         Mxz Mxx Mxy Mx
         Myz Mxy Myy My
         Mz  Mx  My  n]

    F = eigen(M,C)

    values = F.values
    values[values .< 0] .= Inf
    i = argmin(values)

    A,B,C,D = F.vectors[:,i]

    a = -B/(2*A)
    b = -C/(2*A)
    r = sqrt((B^2+C^2-4*A*D)/(4*A^2))
    (a,b,r)
end


end # module
