
# asymmetric test kasa method
for _ in 1:10
    local x0 = 10*rand()
    local y0 = 10*rand()
    local r = rand()+1
    local x = r.*[-1.0,0,0,1,sqrt(2)/2,sqrt(2)/2] .+ x0
    local y = r.*[0.0,1,-1,0,sqrt(2)/2,-sqrt(2)/2] .+ y0
    coefs = CircleFit.kasa(x,y)
    @test x0 ≈ coefs[1]
    @test y0 ≈ coefs[2]
    @test r ≈ coefs[3]
end

# symmetric test kasa method
for _ in 1:10
    local x0 = 10*rand()
    local y0 = 10*rand()
    local r = rand()+1
    local x = r.*[-1.0,0,0,1] .+ x0
    local y = r.*[0.0,1,-1,0] .+ y0
    coefs = CircleFit.kasa(x,y)
    @test x0 ≈ coefs[1]
    @test y0 ≈ coefs[2]
    @test r ≈ coefs[3]
end