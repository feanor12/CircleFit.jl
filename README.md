# CircFit.jl

Circle fitting WIP

Default method used the loss function of ![Kåsa's method](https://doi.org/10.1109/TIM.1976.6312298) combined with lmfit from ![LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) (Levenberg-Marquardt algorithm).

Possible inputs:
* Array of x, Array of y, initial guess
* Array of x, Array of y, Array of weights, initial guess
* Matrix of weights, initial guess

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
