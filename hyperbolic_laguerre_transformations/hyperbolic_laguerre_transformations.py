from PIL import Image, ImageDraw
from numpy import eye, block, sin, cos, tan, arctan2, sign, array, pi, kron, set_printoptions
from numpy.linalg import matrix_power
from scipy.linalg import logm, expm, sqrtm, inv, det
try:
    from .display import display
except ImportError:
    from display import display

set_printoptions(suppress=True)

one = eye(2)
zero = 0*one
jay = array([[0,1],[1,0]])

def double_number(a,b):
    """Uses a 2x2 matrix to represent a double number."""
    return a*one + b*jay

def double_ratio(a,b,c,d):
    """Uses a 4x2 matrix to represent a pair of double numbers, understood
    as a point on the projective double line."""
    return block([[double_number(a,b)],
                  [double_number(c,d)]])

def double_matrix(A,B,C,D,E,F,G,H):
    """Uses a 4x4 matrix to represent a 2x2 matrix of double numbers."""
    return block([
        [A*one+B*jay,C*one+D*jay],
        [E*one+F*jay,G*one+H*jay]])

def make_line(start_theta, end_theta):
    """Generates a double number ratio that represents a line, where
    theta is the angle made with the x-axis, and R is the perpendicular
    distance from the origin."""
    startc = cos(start_theta/2)
    starts = sin(start_theta/2)
    endc = cos(end_theta/2)
    ends = sin(end_theta/2)
    return double_ratio((starts - endc)/2,(starts + endc)/2,
                        (startc + ends)/2,(startc - ends)/2)

def matrix_pair(m1,m2):
    m1_as_double_matrix = kron(m1,eye(2))
    m2_as_double_matrix = kron(m2,eye(2))
    jay_as_double_matrix = double_matrix(0,1,0,0,0,0,0,1)
    return m1_as_double_matrix @ (eye(4) + jay_as_double_matrix)/2 + m2_as_double_matrix.transpose() @ (eye(4) - jay_as_double_matrix)/2

def double_hermitian(m1):
    return matrix_pair(m1,m1)

def inv_sqrt_double_number(double):
    """Takes one over the square root of a double number."""
    return sqrtm(double)

def squared(mat):
    """Multiplies a matrix by itself."""
    return mat @ mat

def normalise(dr):
    """Normalises a double number ratio [a:b] to [a':b'] such that
    (a')**2 + (b')**2 == 1"""
    normalisation_constant = inv_sqrt_double_number(squared(dr[0:2,:]) +
                                                    squared(dr[2:4,:]))
    return dr @ normalisation_constant

def get_line(dr):
    """Returns the (start_theta,end_theta) values of a line 'dr' where
    which give the starting position and ending position of the line
    on the boundary of the hyperbolic plane."""
    # Do the following line twice, because otherwise it somehow fails to
    # normalise it some of the time
    dr = normalise(dr)
    dr = dr.real

    starts = dr[0,0] + dr[0,1]
    endc = dr[0,1] - dr[0,0]
    startc = dr[2,0] + dr[2,1]
    ends = dr[2,0] - dr[2,1]

    return arctan2(starts,startc)*2, arctan2(ends,endc)*2

def double_matrix_determinant(m):
    a = m[0:2,0:2]
    b = m[0:2,2:4]
    c = m[2:4,0:2]
    d = m[2:4,2:4]
    return a @ d  - b @ c

class DeterminantNegativeException(Exception):
    pass

def interpolate(transformation, nframes=50):
    """Returns a list of transformations that interpolates between the
    identity transformation and the specified transformation. If 'nframes'
    equals 1, then simply return the input transformation in a singleton
    list."""
    if nframes == 1:
        return [transformation]
    if double_matrix_determinant(transformation)[0,0] <= 0:
        raise DeterminantNegativeException
    log_transformation = logm(transformation)
    return [expm(log_transformation * i/nframes).real
            for i in range(nframes)]

def apply_transformations(transformations, lines):
    """Applies a list of transformations to a list of lines, generating
    a list of 'frames'. A 'frame' in this case is a list of lines
    each represented as a double number ratio."""
    frames = []
    nframes= len(transformations)
    for i in range(nframes):
        frames.append([])
        transformation = transformations[i]
        for line in lines:
            frames[-1].append(transformation @ line)
    return frames

def draw_frames(frames, imgx=600, imgy=600, width=1, method='poincare', nsamples=50):
    """Returns a list of images which together make up an animation."""
    nframes = len(frames)
    images = [Image.new("RGB", (imgx, imgy)) for i in range(nframes)]
    for i in range(len(frames)):
        draw = ImageDraw.Draw(images[i])
        draw.ellipse((300-100,300-100,300+100,300+100),outline=50,width=3)
        for line in frames[i]:
            try:
                start_theta, end_theta = get_line(line)
            except Exception:
                print(line)
                break
            if method == 'klein':
                draw.line((cos(start_theta)*100+300, sin(start_theta)*100+300,
                        cos(end_theta)*100+300, sin(end_theta)*100+300),
                        width=width, fill=128)
            elif method == 'poincare':
                startx, starty = 300 + 100*cos(start_theta), 300 + 100*sin(start_theta)
                twice_arc = find_centre_of_orthogonal_circle_and_twice_rotation(300,300,300 + 100*cos(start_theta),300 + 100*sin(start_theta),300 + 100*cos(end_theta),300 + 100*sin(end_theta))
                number_of_arcs = 10
                log_overall_arc = number_of_arcs * logm(twice_arc) / 2
                nsamples = 50
                arc_increment = expm(log_overall_arc / nsamples)
                for i in range(nsamples):
                    # this_sample = matrix_power(arc_increment, i)#expm(log_twice_arc * i / nsamples)
                    if i == 0:
                        this_sample = expm(-log_overall_arc / 2)
                    else:
                        this_sample = next_sample
                    next_sample = this_sample @ arc_increment #matrix_power(arc_increment, i+1)#expm(log_overall_arc * (i+1) / nsamples)
                    start_vector = get_coordinates_from_point(this_sample @ point(startx, starty) @ inv(this_sample))
                    end_vector = get_coordinates_from_point(next_sample @ point(startx, starty) @ inv(next_sample))
                    draw.line((start_vector[0],start_vector[1],end_vector[0],end_vector[1]),width=width, fill=128)
    return images

def animate_transformation_without_display(transformation,
                                           lines,
                                           nframes=100,
                                           offset=(0,0),
                                           width=1,
                                           title=None):
    try:
        intermediate_transformations = interpolate(transformation, nframes=nframes)
    except DeterminantNegativeException:
        print("Determinant of the input matrix is negative. Call instead with nframes=1.")
        print("Note that this determinant is double-number valued.")
        return
    frames = apply_transformations(intermediate_transformations, lines)
    images = draw_frames(frames, offset=offset, width=width)
    return images


def animate_transformation(transformation,
                           lines,
                           nframes=100,
                           offset=(0,0),
                           width=1,
                           title=None,
                           via='tk_multiprocess'):
    """Takes as input a transformation and a list of lines.
    It interpolates the transformation starting from the identity
    transformation. It then animates the result of applying the sequence
    of interpolated transformations to each line, and displays the result."""
    try:
        intermediate_transformations = interpolate(transformation, nframes=nframes)
    except DeterminantNegativeException:
        print("Determinant of the input matrix is negative. Call instead with nframes=1.")
        print("Note that this determinant is double-number valued.")
        return
    frames = apply_transformations(intermediate_transformations, lines)
    images = draw_frames(frames, offset=offset, width=width)
    return display(images, title, via=via)


### 2D PGA stuff


q1 = eye(4)
eps = array([[0,1],[0,0]])
epsilon = block([[eps,zero],[zero,eps]])
qi = block([[zero,-one],[one,zero]])
qj = block([[1j*one,zero],[zero,-1j*one]])
qk = qi @ qj

def hodge(q):
    [[ e,  f,  g,  h],
       [   _,  e,    _,  g],
       [ g,  h, _, _],
       [   _,  _,    _, _]] = q.imag
    [[ a,  b, _, _],
       [   _,  a,    _, _],
       [ c,  d,  a,  b],
       [   _,  c,    _,  a]] = q.real
    return a*q1@epsilon + b*q1 + c*qi@epsilon + d*qi + e*qj@epsilon + f*qj + g*qk@epsilon + h*qk

def point(x,y):
    return qi + epsilon @ (x*qj + y*qk)

def cross(q,r):
    return (q@r - r@q)/2

def find_centre_of_orthogonal_circle_and_twice_rotation(x0,y0,x1,y1,x2,y2):
    line_origin_to_1 = hodge(cross(hodge(point(x0,y0)),hodge(point(x1,y1))))
    line_origin_to_2 = hodge(cross(hodge(point(x0,y0)),hodge(point(x2,y2))))
    ninety_from_point_1 = (q1 + point(x1, y1))/sqrt(2)
    ninety_from_point_2 = (q1 + point(x2, y2))/sqrt(2)
    line_1_to_centre = ninety_from_point_1 @ line_origin_to_1 @ inv(ninety_from_point_1)
    line_2_to_centre = ninety_from_point_2 @ line_origin_to_2 @ inv(ninety_from_point_2)
    unnormalised_output = line_2_to_centre @ line_1_to_centre
    return unnormalised_output / det(unnormalised_output)**0.25

def get_coordinates_from_point(q):
    [[ e,  f,  g,  h],
       [   _,  e,    _,  g],
       [ g,  h, _, _],
       [   _,  _,    _, _]] = q.imag
    return (f,h)