

x0 = 1
y0 = 2
n = 20
r0 = 5 .+ randn(n)/10
p = r0 .*exp.(im.*rand(0:0.1:2π,n))
x = real(p) .+ x0 
y = imag(p) .+ y0 
(a,b,r) = CircleFit.pratt(x,y)
@test isapprox(a,x0,atol=0.1)
@test isapprox(b,y0,atol=0.1)
@test isapprox(r,5,atol=0.1)

(a,b,r) = CircleFit.pratt_newton(x,y)
@test isapprox(a,x0,atol=0.1)
@test isapprox(b,y0,atol=0.1)
@test isapprox(r,5,atol=0.1)