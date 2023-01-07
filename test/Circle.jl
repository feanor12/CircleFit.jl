
x0 = 10
y0 = 10
r = 1
x = r.*[-1.1,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0
result = Circle([x0,y0],r,[x y],:kasa)
coefs = coef(result)
@test x0 ≈ coefs[1]
@test y0 ≈ coefs[2]
@test r ≈ coefs[3]
@test dof(result) == 1
@test all(coef(result) .== [result.position...,result.radius])
@test length(coefnames(result)) == length(coef(result))
@test residuals(result) ≈ [0.1,0,0,0]
@test rss(result) ≈ 0.1^2


x = r.*[-1,0,0,1] .+ x0
result = fit(Circle,x,y)
coefs = coef(result)
@test x0 ≈ coefs[1]
@test y0 ≈ coefs[2]
@test r ≈ coefs[3]
@test algorithm(result) == :kasa

result = fit(Circle,x,y,alg=:pratt)
coefs = coef(result)
@test x0 ≈ coefs[1]
@test y0 ≈ coefs[2]
@test r ≈ coefs[3]
@test algorithm(result) == :pratt

dump(x)
dump(y)
result = fit(Circle,x,y,alg=:taubin)
dump(result)
coefs = coef(result)
@test x0 ≈ coefs[1]
@test y0 ≈ coefs[2]
@test r ≈ coefs[3]
@test algorithm(result) == :taubin
