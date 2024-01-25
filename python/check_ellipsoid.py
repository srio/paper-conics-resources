import numpy

#
# get coefficients from table 4
#
def ellipsoid(p=1e10, q=10,theta=3e-3):

    return [1,
            numpy.sin(theta)**2,
            1 - (numpy.sin(theta) * (p-q) / (p+q))**2,
            0,
            -2 * numpy.sin(theta) * numpy.cos(theta) * (q - p) / (p + q),
            0,
            0,
            0,
            -4 * numpy.sin(theta) * p * q / (p + q),
            0,
            ]

#
# get coefficients from table 5
#
def ellipsoid2(input):
    X, N, n, a, b, c = input
    b_over_a = b / a
    b_over_a_square = b_over_a**2
    return [1,
            n[1]**2 + b_over_a_square * n[2]**2,
            n[2]**2 + b_over_a_square * n[1]**2,
            0,
            -2 * n[1] * n[2] * (1 - b_over_a_square),
            0,
            0,
            0,
            2 * (n[2] * X[2] + b_over_a_square * n[1] * X[1]),
            0 ]

def ellipsoid_coordinates(p=1e10, q=10, theta=3e-3): # table 3
    a = (p + q) / 2
    b = numpy.sqrt(p * q) * numpy.sin(theta)
    c = numpy.sqrt(a**2 - b**2)
    Yc = (p**2 - q**2) / (4 * c)
    X = numpy.array([0, Yc, -b * numpy.sqrt(1- (Yc/a)**2)])
    N = numpy.array([0, -2 * Yc / a**2, -2 * X[2] / b**2])
    n = N / numpy.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
    return X, N, n, a, b, c

#
# main
#
def ellipsoid_check(p=200.0, q=10000.0, theta=1e-3, do_assert=1):
    print("\n\n=============================================")
    print("p,q,theta: ",p,q,theta)
    c_p_1 = ellipsoid(p, q, theta)
    c_p_2 = ellipsoid2(ellipsoid_coordinates(p, q, theta))

    for i in range(10):
        diff = numpy.abs(c_p_1[i] - c_p_2[i])
        if  diff < 1e-5:
            ss = ''
        else:
            ss='<<<<<  PROBLEM >>>>>'
        print(i, c_p_1[i], c_p_2[i], ss)
        if do_assert: assert(diff < 1e-5)

    X, N, n, a, b, c = ellipsoid_coordinates(p, q, theta)
    print('X: ', X, [0, (p**2 - q**2)/(4 * c), -b**2/(c*numpy.tan(theta))])
    print('N: ', N, [0, (q**2 - p**2) / (2 * a**2 * c), 2 / (c * numpy.tan(theta))])
    print('n: ', n)

if __name__ == "__main__":
    ellipsoid_check(p=200.0, q=10.0, theta=1e-3, do_assert=1)
    ellipsoid_check(p=10, q=200, theta=1e-3, do_assert=1)

    #
    # random tests
    #
    ntimes = 100
    for i in range(ntimes):
        p = 1000 * numpy.random.rand()
        q = 1000 * numpy.random.rand()
        theta = numpy.random.rand()
        ellipsoid_check(p=p, q=q, theta=theta, do_assert=1)
