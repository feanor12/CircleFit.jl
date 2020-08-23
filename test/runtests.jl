using Test

using CircFit

# symmetric test
for _ in 1:10
    x0 = rand()
    y0 = rand()
    r = rand()
    x = r.*[-1.0,0,0,1] .+ x0
    y = r.*[0.0,1,-1,0] .+ y0
    p0 = [0.0,0,1]
    result = circfit(x,y,p0)
    @assert x0 ≈ result.param[1]
    @assert y0 ≈ result.param[2]
    @assert r ≈ result.param[3]
end

# asymmetric test
for _ in 1:10
    x0 = rand()
    y0 = rand()
    r = rand()
    x = r.*[-1.0,0,0,1,sqrt(2)/2,sqrt(2)/2] .+ x0
    y = r.*[0.0,1,-1,0,sqrt(2)/2,-sqrt(2)/2] .+ y0
    p0 = [0.0,0,1]
    result = circfit(x,y,p0)
    @assert x0 ≈ result.param[1]
    @assert y0 ≈ result.param[2]
    @assert r ≈ result.param[3]
end
