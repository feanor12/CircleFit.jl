# Provided by Dennis Hoelgaard Bal KronosTheLate 
# https://github.com/feanor12/CircleFit.jl/issues/10

using LinearAlgebra

# Simplified helper function to extract parameters from the multivector
function extract_mv_parameters2(mv::Vector{<:Real})
    cvx, cvy, cvz = mv[1:3]   # e2^e3, e3^e1, e1^e2
    cgx, cgy, cgz = mv[4:6]   # e0^e1, e0^e2, e0^e3
    cmx, cmy, cmz = mv[7:9]   # e1^einf, e2^einf, e3^einf
    cgw = -mv[10]             # e0^einf

    num = cvx^2 + cvy^2 + cvz^2 - cgw^2 + 2 * (cgx * cmx + cgy * cmy + cgz * cmz)
    denom = cgx^2 + cgy^2 + cgz^2
    r = sqrt(num / denom)

    n = -[cgx, cgy, cgz]
    nn = n / norm(n)

    B = [cgw -cvz cvy; cvz cgw -cvx; -cvy cvx cgw]
    c = (B * n) / sum(n .^ 2)

    return c, nn, r
end

# Main function to fit a circle in 3D space using CGA
"""
Takes in points in 3D space and outputs a fitted circle in 3D based on a method in Conformal Geometric Algebra (CGA).

    Input is points in 3D euclidean space with each row being a point, and each column is the xyz coordinate.
    The ouput is: centrum, radius, normal_vector
    where circle_cga is the CGA representation of the fitted circle, in the IPNS representation!

    This method is solved by calculating the stacked matrix: M_stacked, which is solved with SVD.

    Example usage:
    >>> centrum, radius, normal_vector, circle_cga = circlefitCGA3D(points)
"""
function circlefitCGA3D(points)
    N = size(points, 1)
    M_stacked = zeros(Float64, 5 * N, 10)

    for i in 1:N
        x4 = 0.5 * sum(points[i, :] .^ 2)

        px, py, pz = points[i, :]
        vec_cross = [0 -pz py; pz 0 -px; -py px 0]

        Mi = [
            -vec_cross  -x4 * I(3)  I(3)  zeros(3, 1);
            zeros(1, 3)  -points[i, :]'  zeros(1, 3)  1;
            zeros(1, 3)  zeros(1, 3)  points[i, :]'  -x4
        ]

        M_stacked[5*(i-1)+1:5*i, :] = Mi
    end

    V = svd(M_stacked, full=false).V
    circle_og = V[:, end]

    return extract_mv_parameters2(circle_og)
end

# Wrapper function to accept separate x, y, z vectors
circlefitCGA3D(xs, ys, zs) = circlefitCGA3D(hcat(xs, ys, zs))

