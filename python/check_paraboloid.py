import numpy

#
# get coefficients from table 4
#
def paraboloid(p=1e10, q=10,theta=3e-3):
    if p > q:
        return paraboloid_focusing(q=q, theta=theta)
    else:
        return paraboloid_collimating(p=p, theta=theta)

def paraboloid_focusing(q=10, theta=3e-3):
    print('focusing')
    return [1,
            numpy.sin(theta)**2,
            numpy.cos(theta)**2,
            0,
            2 * numpy.cos(theta) * numpy.sin(theta),
            0,
            0,
            0,
            -4 * q * numpy.sin(theta),
            0,
            ]

# see conics_penelope_paraboloid_collimating.nb
def paraboloid_collimating(p=10, theta=3e-3):
    print('collimating')
    return [1,
            numpy.sin(theta)**2,
            numpy.cos(theta)**2,
            0,
            -2 * numpy.cos(theta) * numpy.sin(theta),
            0,
            0,
            0,
            -4 * p * numpy.sin(theta),
            0,
            ]

#
# get coefficients from table 5
#
def paraboloid2(input):
    X, N, n, a_p = input
    return [1, n[2]**2, n[1]**2, 0, 2 * n[1] * n[2], 0, 0, 0, -4 * a_p / n[2], 0 ]

def paraboloid_coordinates(p=1e10, q=10,theta=3e-3): # table 3
    if p > q:
        return paraboloid_coordinates_focusing(q=q, theta=theta)
    else:
        return paraboloid_coordinates_collimating(p=p, theta=theta)

def paraboloid_coordinates_focusing(q=10, theta=3e-3):
    print('focusing')
    a_p = q * numpy.sin(theta)**2
    X = numpy.array([0, -q * numpy.sin(2 * theta), q * numpy.cos(theta)**2])
    N = numpy.array([0, 2 * q * numpy.sin(2 *theta), 4 * a_p])
    n = N / numpy.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
    return X, N, n, a_p

def paraboloid_coordinates_collimating(p=10, theta=3e-3):
    print('collimating')
    a_p = p * numpy.sin(theta)**2
    X = numpy.array([0, p * numpy.sin(2 * theta), p * numpy.cos(theta)**2])
    N = numpy.array([0, -2 * p * numpy.sin(2 * theta), 4 * a_p])
    n = N / numpy.sqrt(N[0]**2 + N[1]**2 + N[2]**2)
    return X, N, n, a_p


#
# main
#
def paraboloid_check(p=200.0, q=10000.0, theta=1e-3, do_assert=1):
    print("\n\n=============================================")
    print("p,q,theta: ",p,q,theta)
    c_p_1 = paraboloid(p, q, theta)
    c_p_2 = paraboloid2(paraboloid_coordinates(p, q, theta))

    for i in range(10):
        diff = numpy.abs(c_p_1[i] - c_p_2[i])
        if  diff < 1e-5:
            ss = ''
        else:
            ss='<<<<<  PROBLEM >>>>>'
        print(i, c_p_1[i], c_p_2[i], ss)
        if do_assert: assert(diff < 1e-5)

    X, N, n, a_p = paraboloid_coordinates(p, q, theta)
    if p > q: # focusing
        print('X: ', X, [0, -2 * q * numpy.sin(theta) * numpy.cos(theta), q * numpy.cos(theta) ** 2])
    else: # colllimatng
        print('X: ', X, [0, 2 * p * numpy.sin(theta) * numpy.cos(theta), p * numpy.cos(theta) ** 2])
    print('N: ', N, [0, -2 * X[1], 4 * a_p])
    print('n: ', n, [0, numpy.cos(theta), numpy.sin(theta)])

if __name__ == "__main__":
    paraboloid_check(p=200.0, q=10.0, theta=1e-3, do_assert=1)
    paraboloid_check(p=10, q=200, theta=1e-3, do_assert=1)

    #
    # random tests
    #
    ntimes = 100
    for i in range(ntimes):
        p = 1000 * numpy.random.rand()
        q = 1000 * numpy.random.rand()
        theta = numpy.random.rand()
        paraboloid_check(p=p, q=q, theta=theta, do_assert=1)
