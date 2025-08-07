using LinearAlgebra

x0 = 1
y0 = 2
n = 20
r0 = 5 .+ randn(n)/10
p = r0 .*exp.(im.*rand(0:0.1:2π,n))
x = real(p) .+ x0 
y = imag(p) .+ y0 
c =  StatsBase.fit(Circle,x,y,zeros(size(x)))
@test isapprox(c.position[1],x0,atol=0.1)
@test isapprox(c.position[2],y0,atol=0.1)
@test isapprox(c.radius,5,atol=0.1)


"""
Generates points on a circle in 3D and outputs the 3D euclidean coordinates of the generated points.

The function takes in a normal vector, which expresses the direction of the circle, meaning the direction of the plane the circle lies in,
    radius of the circle, offsets to displace the circle center in 3D, how many points are to be generated, and how large the random noise to be added is.
    The output is a matrix of size num_points x 3, with the format:
    [[x1, y1, z1],
    [x2, y2, z2],
    [x3, y3, z3]...]

    The noise is independent Gaussian noise on the x, y, and z-coordinate.

    Example usage:
    >>> circle_points = generate_circle_points([1.0, 1.0, 1.0], 2.0, [1.0, 1.0, 1.0], 100, 0.05)
"""
    function generate_circle_points(normal_vector::Vector{<:Real}, radius::Real, offsets::Vector{<:Real}, num_points::Integer, noise_scale::Real)
    # Generate points on a 2D circle
    theta = range(0, stop=2π*(num_points-1)/num_points, length=num_points)
    x_circle = radius * cos.(theta)
    y_circle = radius * sin.(theta)

    # Normalize the normal vector
    normal_vector /= norm(normal_vector)

    # Calculate the axis-angle representation for rotation
    # Reference: https://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
    angle = acos(normal_vector[3])  # Angle between z-axis and normal vector
    axis = cross([0.0, 0.0, 1.0], normal_vector)  # Axis of rotation
    axis /= norm(axis)

    # Rotation matrix
    c = cos(angle)
    s = sin(angle)
    t = 1 - c
    rotation_matrix = [
        t*axis[1]^2 + c                 t*axis[1]*axis[2] - s*axis[3]   t*axis[1]*axis[3] + s*axis[2]
        t*axis[1]*axis[2] + s*axis[3]   t*axis[2]^2 + c                 t*axis[2]*axis[3] - s*axis[1]
        t*axis[1]*axis[3] - s*axis[2]   t*axis[2]*axis[3] + s*axis[1]   t*axis[3]^2 + c
    ]

    # Rotate points to align with the specified normal vector
    rotated_points = rotation_matrix * hcat(x_circle, y_circle, zeros(num_points))'
    x_circle, y_circle, z_circle = rotated_points[1, :], rotated_points[2, :], rotated_points[3, :]

    # Set random seed for reproducibility of 'random' noise. Can be commented out.
    # Random.seed!(42)

    # Add noise to the points and the constant offset
    x_circle += randn(num_points) * noise_scale .+ offsets[1]
    y_circle += randn(num_points) * noise_scale .+ offsets[2]
    z_circle += randn(num_points) * noise_scale .+ offsets[3]

    return hcat(x_circle, y_circle, z_circle)
end


let
    # Example usage
    normal_vector = [1, 1, 1]  # Normal vector
    radius_og = 1
    num_points = 100
    offsets = [1, -1, 0.25]

    noise_scale = 0.05
    points = generate_circle_points(normal_vector, radius_og, offsets, num_points, noise_scale)


    # 'Default method for circle fitting, which is SVD on M_stacked
    # centrum, radius, normal_vector = circlefitCGA3D(points)
    centrum, normal_vector_out, radius = CircleFit.circlefitCGA3D(points)
    @test all(isapprox.(centrum,offsets,atol=0.1))
    # proj to normal vector
    @test isapprox(abs(normal_vector_out'*(normal_vector./sqrt(3))),1,atol=0.1)
    @test isapprox(radius,radius_og,atol=0.1)
end
