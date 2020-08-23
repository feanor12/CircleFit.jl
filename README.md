# Circfit

Circle fitting WIP

Default method used the loss function of Ka˙sa's method combined with lmfit from LsqFit (Levenberg-Marquardt algorithm).

Example:

```julia
using Circfit
r = 5
x0 = 2
y0 = 4.5
x = r.*[-1.0,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0
p0 = [0.0,0,1]
result = circfit(x,y,p0)
println("x0,y0,r : ", result.param)
# x0,y0,r : [1.9999999999999991, 4.499999999999997, 5.00000000000002]
```