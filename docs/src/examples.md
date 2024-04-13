
```@eval
using StatsBase
using CircleFit
using Plots

for alg in [:kasa,:pratt,:taubin,:graf]
    plots = map(Iterators.product(0.5:0.5:2,[5,10,20])) do (phi_max,npoints)
        actual_points = round(Int,npoints*phi_max)
        t = range(0,pi*phi_max-2*pi/(actual_points),length=actual_points)
        x = sin.(t) .+ 2 .+ rand(size(t)...)/10
        y = cos.(t) .+ 1 .+ rand(size(t)...)/10

        p=scatter(x, y, color="red",aspect_ratio=1,
            ylims=(-0.5,2.5),xlims=(0.5,3.5),legend=false)
        circle = fit(Circle,x,y,alg=alg);
        x_fit,y_fit = CircleFit.parametric_form(circle);
        plot!(x_fit,y_fit,label=string(alg),color=:black);

        return p
    end
    p = plot(plots...,plot_title=string(alg),layout=(3,4))
    savefig(p,"plot_$(string(alg)).svg")
end
nothing
```
![](plot_kasa.svg)
![](plot_pratt.svg)
![](plot_taubin.svg)
![](plot_graf.svg)
