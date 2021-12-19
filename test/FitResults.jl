
x0 = 10
y0 = 10
r = 1
x = r.*[-1.1,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0
result = CircleFit.FitResult([x0,y0],r,[x y])
coefs = coef(result)
@test x0 ≈ coefs[1]
@test y0 ≈ coefs[2]
@test r ≈ coefs[3]
@test dof(result) == 1
@test all(coef(result) .== [result.position...,result.radius])
@test length(coefnames(result)) == length(coef(result))
@test residuals(result) ≈ [0.1,0,0,0]
@test rss(result) ≈ 0.1^2
