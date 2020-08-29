
# asymmetric test kasa method
for _ in 1:10
    x0 = 10*rand()
    y0 = 10*rand()
    r = rand()+1
    x = r.*[-1.0,0,0,1,sqrt(2)/2,sqrt(2)/2] .+ x0
    y = r.*[0.0,1,-1,0,sqrt(2)/2,-sqrt(2)/2] .+ y0
    result = circfit(x,y)
    @assert x0 ≈ result[1]
    @assert y0 ≈ result[2]
    @assert r ≈ result[3]
end

# symmetric test kasa method
for _ in 1:10
    x0 = 10*rand()
    y0 = 10*rand()
    r = rand()+1
    x = r.*[-1.0,0,0,1] .+ x0
    y = r.*[0.0,1,-1,0] .+ y0
    result = circfit(x,y)
    @assert x0 ≈ result[1]
    @assert y0 ≈ result[2]
    @assert r ≈ result[3]
end