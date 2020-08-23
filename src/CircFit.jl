module CircFit
import LsqFit: lmfit
import Statistics: var, cov

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
    x² = x.^2
    y² = y.^2
    
    A = n*(n-1)*var(x) 
    B = n*(n-1)*cov(x,y) 
    C = n*(n-1)*var(y) 
    D = n*(n-1)/2*(cov(x,y²)+cov(x,x²)) 
    E = n*(n-1)/2*(cov(y,x²)+cov(y,y²)) 

    ACB2 = (A*C-B^2)
    am = (D*C - B*E) / ACB2
    bm = (A*E - B*D) / ACB2
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
