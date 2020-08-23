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


# test weights
for _ in 1:10
    x0 = rand()
    y0 = rand()
    r = rand()
    x = r.*[-1.0,0,0,1,sqrt(2)/2,sqrt(2)/2,0,0,0] .+ x0
    y = r.*[0.0,1,-1,0,sqrt(2)/2,-sqrt(2)/2,1.1,1.1,1.1] .+ y0
    p0 = [0.0,0,1]
    result1 = circfit(x,y,p0)
    x = r.*[-1.0,0,0,1,sqrt(2)/2,sqrt(2)/2,0] .+ x0
    y = r.*[0.0,1,-1,0,sqrt(2)/2,-sqrt(2)/2,1.1] .+ y0
    wt = [1,1,1,1,1,1,3]
    result2 = circfit(x,y,wt,p0)
    wt = wt ./ sum(wt)
    result3 = circfit(x,y,wt,p0)

    for i in 1:3
        @assert result1.param[i] ≈ result2.param[i] ≈ result3.param[i]
    end
end
