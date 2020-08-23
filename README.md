# Circfit

Circle fitting WIP

Default method used the loss function of Ka˙sa's method combined with lmfit from LsqFit (Levenberg-Marquardt algorithm).

```math
loss(xdata,ydata,x0,y0,r) = \sum_i \left( r^2 - (x0 - xdata_i)^2 -(y0-ydata_i)^2 \right)^2
```

Example:

```julia
# import library
using Circfit
# generate test data
r = 5
x0 = 2
y0 = 4.5
x = r.*[-1.0,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0
# initial fit parameters
p0 = [0.0,0,1]
# fit
result = circfit(x,y,p0)
# print results
println("x0,y0,r : ", result.param)
# x0,y0,r : [1.9999999999999991, 4.499999999999997, 5.00000000000002]
```
