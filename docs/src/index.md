# CircleFit Documentation

This package includes algorithms used for fitting circles.

| algorithm | 2D | 3D |source|
|-----------|----|----|------|
| kasa| yes | no |[A circle fitting procedure and its error analysis](https://doi.org/10.1109/TIM.1976.6312298) |
| graf| yes | no |[Least Squares Fitting of Circles](https://link.springer.com/article/10.1007/s10851-005-0482-8)|
| taubin| yes | no |[Least Squares Fitting of Circles](https://link.springer.com/article/10.1007/s10851-005-0482-8)|
| pratt| yes | no | [Least Squares Fitting of Circles](https://link.springer.com/article/10.1007/s10851-005-0482-8)|
| pratt\_newton| yes | no | [Circular and Linear Regression](https://doi.org/10.1201/EBK1439835906), [C++ code](https://people.cas.uab.edu/~mosya/cl/CircleFitByPratt.cpp), [Least Squares Fitting of Circles](https://link.springer.com/article/10.1007/s10851-005-0482-8)|

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
