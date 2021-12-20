
@test all(CircleFit.BCD_to_abr(CircleFit.abr_to_BCD(3.0,4.0,6.0)...) .≈ (3.0,4.0,6.0))

x0 = 1.0
y0 = 2.0
n = 20
r0 = 5
p = r0 .*exp.(im.*range(0,2π,length=n))
x = real(p) .+ x0 
y = imag(p) .+ y0 
(a,b,r) = CircleFit.GRAF(x,y,[1.0,1.0,1.0])
@test a ≈ x0
@test b ≈ y0
@test r ≈ r0

(a,b,r) = coef(fit(Circle,x,y,alg=:graf))
@test a ≈ x0
@test b ≈ y0
@test r ≈ r0