import numpy

#
# get coefficients from table 4
#
def hyperboloid(p=1e10, q=10,theta=3e-3):

    return [1,
            numpy.sin(theta)**2,
            1 - (numpy.sin(theta) * (p + q) / (p - q))**2,
            0,
            -2 * numpy.sin(theta) * numpy.cos(theta) * (q + p) / (q - p),
            0,
            0,
            0,
            -4 * numpy.sin(theta) * p * q / (q - p),
            0,
            ]

#
# get coefficients from table 5
#
def hyperboloid2(input):
    X, N, n, a, b, c = input
    b_over_a = b / a
    b_over_a_square = b_over_a**2

    # A = -1/b**2
    # B = 1/a**2
    # out = [A,
    #         A * n[1]**2 + B * n[2]**2,
    #         A * n[2]**2 + B * n[1]**2,
    #         0,
    #         2 * n[1] * n[2] * (B - A),
    #         0,
    #         0,
    #         0,
    #         2 * (A * n[2] * X[2] - B * n[1] * X[1]),
    #         0 ]
    #
    # out = numpy.array(out)
    # return out / A

    return [1,
            n[1]**2 - b_over_a_square * n[2]**2,
            n[2]**2 - b_over_a_square * n[1]**2,
            0,
            -2 * n[1] * n[2] * (1 + b_over_a_square),
            0,
            0,
            0,
            2 * (n[2] * X[2] - b_over_a_square * n[1] * X[1]),
            0 ]

def ellipsoid_coordinates(p=1e10, q=10, theta=3e-3): # table 3
    a = numpy.abs(p - q) / 2
    c = 0.5 * numpy.sqrt(p**2 + q**2 - 2 * p * q * numpy.cos(2 * theta))
    b = numpy.sqrt(c**2 - a**2)
    Yc = (p**2 - q**2) / (4 * c)
    X = numpy.array([0, Yc, b * numpy.sqrt((Yc/a)**2 - 1)])
    N = numpy.array([0, -2 * Yc / a**2, 2 * X[2] / b**2])
    if q > p:
        N *= -1
    n = N / numpy.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
    return X, N, n, a, b, c

#
# main
#
def hyperboloid_check(p=200.0, q=10000.0, theta=1e-3, do_assert=1):
    print("\n\n=============================================")
    print("p,q,theta: ",p,q,theta)
    c_p_1 = hyperboloid(p, q, theta)
    c_p_2 = hyperboloid2(ellipsoid_coordinates(p, q, theta))

    for i in range(10):
        diff = numpy.abs(c_p_1[i] - c_p_2[i])
        if  diff < 1e-5:
            ss = ''
        else:
            ss='<<<<<  PROBLEM >>>>>'
        print(i, c_p_1[i], c_p_2[i], ss)
        if do_assert: assert(diff < 1e-5)

    X, N, n, a, b, c = ellipsoid_coordinates(p, q, theta)
    print('X: ', X, [0, (p**2 - q**2)/(4 * c), b * numpy.sqrt(X[1]**2/a**2 - 1)])
    if p > q:
        print('N: ', N, [0, -2 * X[1] / a**2,  2 * X[2] / b**2])
    else:
        print('N: ', N, [0,  2 * X[1] / a**2, -2 * X[2] / b**2])

    print('n: ', n)
    return c_p_1, c_p_2

def height(ccc,y=0,x=0,return_solution=0):

    aa = ccc[2]
    bb = ccc[4] * y + ccc[5] * x + ccc[8]
    cc = ccc[0] * x**2 + ccc[1] * y**2 + ccc[3] * x * y + \
        ccc[6] * x + ccc[7] * y + ccc[9]

    if aa != 0:
        discr = bb**2 - 4 * aa * cc + 0j
        s1 = (-bb + numpy.sqrt(discr)) / 2 / aa
        s2 = (-bb - numpy.sqrt(discr)) / 2 / aa

        if return_solution == 0: # select the solution close to zero at pole
            if numpy.abs(s1).min() < numpy.abs(s2).min():
                ss = s1
            else:
                ss = s2
        elif return_solution == 1:
            ss = s1
        else:
            ss = s2
    else:
        ss = -cc / bb

    return numpy.real(ss)

if __name__ == "__main__":
    # hyperboloid_check(p=200.0, q=10.0, theta=1e-3, do_assert=0)

    p = 200
    q = 10

    p = 10
    q = 200
    c_p_1, c_p_2 = hyperboloid_check(p=p, q=q, theta=1e-3, do_assert=1)



    y = numpy.linspace(-0.05, 0.05, 100)
    z1 = height(c_p_1, y)
    z2 = height(c_p_2, y)

    from srxraylib.plot.gol import plot
    plot(y, z1,
         y, z2,
         legend=['ken','manolo'], title="p: %f, q: %f" % (p,q))


    #
    # random tests
    #
    ntimes = 100
    for i in range(ntimes):
        p = 100 * numpy.random.rand()
        q = 100 * numpy.random.rand()
        theta = numpy.random.rand()
        print("\n\n===================================================", p, q)
        hyperboloid_check(p=p, q=q, theta=theta, do_assert=1)
