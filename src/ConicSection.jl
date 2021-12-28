abstract type ConicSection end

struct Circle <: ConicSection
    radius
    center
end

struct Ellipse <: ConicSection
end

struct Parabola <: ConicSection
end

struct Hyperbola <: ConicSection
end

