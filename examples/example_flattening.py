from laguerre_transformations import *

my_lines = [make_line(t, t + pi/2) for t in linspace(0,2*pi,100)]

transformation = double_hermitian(array([[20,0],[0,1]]))

animate_transformation(transformation,
                       my_lines,
                       offset=(200,200))