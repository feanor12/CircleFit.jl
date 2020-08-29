# CircleFit.jl

[![Build Status](https://travis-ci.com/feanor/CircleFit.jl.svg?branch=master)](https://travis-ci.com/feanor/CircleFit.jl)
[![Coverage](https://codecov.io/gh/feanor/CircleFit.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/feanor/CircleFit.jl)


Circle fitting using [Kåsa's method](https://doi.org/10.1109/TIM.1976.6312298)

Example:
```julia
# import library
using CircleFit
# generate test data
r = 5
x0 = 2
y0 = 4.5
x = r.*[-1.0,0,0,1] .+ x0
y = r.*[0.0,1,-1,0] .+ y0
# fit
x0,y0,radius = circfit(x,y)
#(2.0, 4.499999999999999, 5.0)
```

Non optimized implementations:
* Method by Taubin `CircleFit.taubin(x,y)`
* Method by Pratt `CircleFit.pratt(x,y)`
