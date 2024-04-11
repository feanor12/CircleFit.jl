# CircleFit.jl

[![CI](https://github.com/feanor12/CircleFit.jl/actions/workflows/test.yml/badge.svg)](https://github.com/feanor12/CircleFit.jl/actions/workflows/test.yml)
[![Coverage](https://codecov.io/gh/feanor12/CircleFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/feanor12/CircleFit.jl)
[![][docs-development-img]][docs-development-url]

Circle fitting using [Kåsa's method](https://doi.org/10.1109/TIM.1976.6312298)

Example:
```julia
# import library
using CircleFit
using StatsBase

# generate test data
r = 5
x0 = 2
y0 = 4.5
x = r.*[-1.0,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0

# fit
result = fit(Circle,x,y)

# fitted coefficients
coef(result)
# (2.0, 4.499999999999999, 5.0)

# names of the coefficients
coefnames(result)
# ("center position x1", "center position x2", "radius")

# calculate the residuals square sum
rss(result)
# 1.5777218104420236e-30
```

Gradient weighted algebraic fit for a circle using `LsqFit.levenberg_marquardt`: `CircleFit.GRAF(x,y,p0)` or `fit(Circle,x,y,alg=:graf)`

Non optimized implementations:
* Method by Taubin `CircleFit.taubin(x,y)` or `fit(Circle,x,y,alg=:taubin)`
* Method by Pratt `CircleFit.pratt(x,y)` or `fit(Circle,x,y,alg=:pratt)`
