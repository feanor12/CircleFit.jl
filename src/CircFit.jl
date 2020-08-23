module CircFit
import LsqFit: lmfit

export circfit, Kasa

@enum Method begin
    Kasa
end

"""
Using analytical solution of Kasa
https://ieeexplore.ieee.org/document/1246564
x y values are arrays
"""
function circfit(::Type{Val{Kasa}},x::AbstractArray,y::AbstractArray)
    n = length(x)
    #TODO: use std and cov from Statistics
    A = n*sum(x.^2) - (sum(x))^2
    B = n*sum(x.*y) - sum(x)*sum(y)
    C = n*sum(y.^2) - sum(y)^2
    D = (n*sum(x.*y.^2) - sum(x)*sum(y.^2) + n*sum(x.^3) - sum(x)*sum(x.^2))/2
    E = (n*sum(y.*x.^2) - sum(y)*sum(x.^2) + n*sum(y.^3) - sum(y)*sum(y.^2))/2

    am = (D*C - B*E) / (A*C-B^2)
    bm = (A*E - B*D) / (A*C-B^2)  
    rk = sqrt(sum(((x.-am).^2 + (y.-bm).^2)/n))

    (am,bm,rk)
end

"""
x: Array of x-coordinates
y: Array of y-coordinates
p0: Array of initial values [center_x, center_y, radius] 
"""
function circfit(x::AbstractArray,y::AbstractArray,p0::AbstractArray,kwargs...)
    T = eltype(y)
    f = (p) -> @. p[3]^2 - (x-p[1])^2 - (y-p[2])^2
    lmfit(f,p0,T[];kwargs...)
end

"""
x: Array of x-coordinates
y: Array of y-coordinates
w: Array of weights
p0: Array of initial values [center_x, center_y, radius] 
"""
function circfit(x::AbstractArray,y::AbstractArray,wt::AbstractArray,p0::AbstractArray,kwargs...)
    u = sqrt.(wt)
    f = (p) -> @. u*(p[3]^2 - (x-p[1])^2 - (y-p[2])^2)
    lmfit(f,p0,wt;kwargs...)
end

"""
wt: Array of weights
p0: Array of initial values [center_x, center_y, radius] 

returns the center and radius is terms of the array index (starting at one)
center_x: dim=1
center_y: dims=2
"""
function circfit(wt::AbstractMatrix,p0::AbstractArray,kwargs...)
    u = sqrt.(wt)
    x = repeat(1:size(wt,1),1,size(wt,2))
    y = repeat((1:size(wt,2))',size(wt,1),1)
    
    f = (p) -> view((@. u*(p[3]^2 - (x-p[1])^2 - (y-p[2])^2)),:)
    lmfit(f,p0,view(wt,:);kwargs...)
end

end # module
