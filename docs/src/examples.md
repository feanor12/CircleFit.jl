
```@eval
using StatsBase
using CircleFit
using Plots

for alg in [:kasa,:pratt,:taubin,:graf,:cga]
    plots = map(Iterators.product([0.5,1,2],[10,20])) do (phi_max,npoints)
        actual_points = round(Int,npoints*phi_max)
        t = range(0,pi*phi_max,length=actual_points)
        x = sin.(t) .+ 2 .+ rand(size(t)...)/10
        y = cos.(t) .+ 1 .+ rand(size(t)...)/10

        p=scatter(x, y, color="red",aspect_ratio=1,
            ylims=(-0.5,2.5),xlims=(0.5,3.5),legend=false)
        circle = fit(Circle,x,y,alg=alg);
        x_fit,y_fit = CircleFit.parametric_form(circle);
        plot!(x_fit,y_fit,label=string(alg),color=:black);

        return p
    end
    p = plot(plots...,plot_title=string(alg),layout=(2,3))
    savefig(p,"plot_$(string(alg)).svg")
end
nothing
```
![](plot_kasa.svg)
![](plot_pratt.svg)
![](plot_taubin.svg)
![](plot_graf.svg)
![](plot_cga.svg)


```@eval
using StatsBase
using CircleFit
using LinearAlgebra
using Plots

    plots = map(Iterators.product([0.5,1,2],[10,20])) do (phi_max,npoints)
        actual_points = round(Int,npoints*phi_max)
        t = range(0,pi*phi_max,length=actual_points)
        x = sin.(t) .+ 2 .+ rand(size(t)...)/10
        y = cos.(t) .+ 1 .+ rand(size(t)...)/10
        z = zeros(size(x))

        tilt = pi/4
        x = cos.(tilt).*(x.-2).+2 
        z = -sin.(tilt).*(x.-2)

        p=scatter(x, y, z, color="red",aspect_ratio=1,
            ylims=(-0.5,2.5),xlims=(0.5,3.5),zlims=(-1.5,1.5),legend=false)
        circle = fit(Circle,x,y,z,alg=:cga);

        # calc circle plane
        u = [1,0,0]
        u = u - (u'*circle.normal)*circle.normal
        u = u ./ norm(u)
        v = cross(circle.normal,u)

        t = 0:0.1:(2*pi)

        x_fit = circle.radius .* (sin.(t)*u[1] + cos.(t)*v[1]) .+ circle.position[1]
        y_fit = circle.radius .* (sin.(t)*u[2] + cos.(t)*v[2]) .+ circle.position[2]
        z_fit = circle.radius .* (sin.(t)*u[3] + cos.(t)*v[3]) .+ circle.position[3]
        plot!(x_fit,y_fit,z_fit,color=:black);

        return p
    end
    p = plot(plots...,plot_title="cga",layout=(2,3))
    savefig(p,"3dplot_cga.svg")

nothing
```
![](3dplot_cga.svg)
