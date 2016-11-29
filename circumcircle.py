# code from http://codepad.org/KYyfy9CH
from math import cos
from math import sqrt
# Why math? Just 'cos

# Distance between points A and B
def distance(A, B):
    n = len(A)
    assert len(B) == n
    return sum((A[i] - B[i]) ** 2 for i in range(n)) ** 0.5


# Cosine of angle ABC
def cosine(A, B, C):
    a, b, c = distance(B, C), distance(A, C), distance(A, B)
    return (a * a + c * c - b * b) / (2 * a * c)


# Cartesian coordinates of the point whose barycentric coordinates
# with respect to the triangle ABC are [p,q,r]
def barycentric(A, B, C, p, q, r):
    n = len(A)
    assert len(B) == len(C) == n
    s = p + q + r
    p, q, r = p / s, q / s, r / s
    return tuple([p * A[i] + q * B[i] + r * C[i] for i in range(n)])


# Cartesian coordinates of the point whose trilinear coordinates
# with respect to the triangle ABC are [alpha,beta,gamma]
def trilinear(A, B, C, alpha, beta, gamma):
    a = distance(B, C)
    b = distance(A, C)
    c = distance(A, B)
    return barycentric(A, B, C, a * alpha, b * beta, c * gamma)


# Cartesian coordinates of the circumcenter of triangle ABC
def circumcenter(A, B, C):
    cosA = cosine(C, A, B)
    cosB = cosine(A, B, C)
    cosC = cosine(B, C, A)
    return trilinear(A, B, C, cosA, cosB, cosC)


def radius(A, B, C):
    a = distance(B, C)
    b = distance(A, C)
    c = distance(A, B)
    return ((a * b * c)/sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c)))

# circumcenter((1, 2), (3, 3), (6, 0))
