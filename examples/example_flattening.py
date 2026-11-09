from numpy import linspace, pi
from hyperbolic_laguerre_transformations import (
    animate_transformation,
    double_matrix,
    make_line,
)

my_lines = [make_line(t, t + pi/2) for t in linspace(0, 2*pi, 100)]
transformation = double_matrix(20, 0, 0, 0, 0, 0, 1, 0)

animate_transformation(transformation, my_lines)

