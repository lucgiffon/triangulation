# code from http://codepad.org/KYyfy9CH
import math
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import random
# Why math? Just 'cos


class CircumCircle:
    def __init__(self, A, B, C):
        self.center, self.radius = self.draw(A, B, C)
        # self.center = self.circumcenter(A, B, C)
        # self.radius = self.circumradius(A, B, C)

    # Distance between points A and B
    def distance(self, p1, p2):
        n = len(p1)
        assert len(p2) == n
        return sum((p1[i] - p2[i]) ** 2 for i in range(n)) ** 0.5

    # Cosine of angle p1-p2-p3
    def cosine(self, p1, p2, p3):

        d1, d2, d3 = self.distance(p2, p3), self.distance(p1, p3), self.distance(p1, p2)
        return (d1 * d1 + d3 * d3 - d2 * d2) / (2 * d1 * d3)

    # Cartesian coordinates of the point whose barycentric coordinates
    # with respect to the triangle ABC are [p,q,r]
    def barycentric(self, A, B, C, p, q, r):
        n = len(A)
        assert len(B) == len(C) == n
        s = p + q + r
        p, q, r = p / s, q / s, r / s
        return tuple([p * A[i] + q * B[i] + r * C[i] for i in range(n)])


    # Cartesian coordinates of the point whose trilinear coordinates
    # with respect to the triangle ABC are [alpha,beta,gamma]
    def trilinear(self, A, B, C, alpha, beta, gamma):
        a = self.distance(B, C)
        b = self.distance(A, C)
        c = self.distance(A, B)
        return self.barycentric(A, B, C, a * alpha, b * beta, c * gamma)


    # Cartesian coordinates of the circumcenter of triangle ABC
    def circumcenter(self, A, B, C):
        cosA = self.cosine(C, A, B)
        cosB = self.cosine(A, B, C)
        cosC = self.cosine(B, C, A)
        return self.trilinear(A, B, C, cosA, cosB, cosC)


    def circumradius(self, A, B, C):
        a = self.distance(B, C)
        b = self.distance(A, C)
        c = self.distance(A, B)
        return (a * b * c)/math.sqrt((a+b+c) * (b+c-a) * (c+a-b) * (a+b-c))

    def draw(self, A, B, C):
        pts = np.asarray([np.asarray(p) for p in [A, B, C]])
        pts2 = np.dot(pts, pts.T)
        A = np.bmat([[2 * pts2, [[1],
                                 [1],
                                 [1]]],
                     [[[1, 1, 1, 0]]]])

        b = np.hstack((np.sum(pts * pts, axis=1), [1]))

        x = np.linalg.solve(A, b)
        bary_coords = x[:-1]
        center = np.dot(bary_coords, pts)

        radius = np.linalg.norm(pts[0] - center) # euclidean distance
        # radius = math.sqrt(np.sum(np.square(pts[0] - center)))  # squared distance
        return center, radius

    def contains(self, p):
        # using the pythagore theorem, we test if p lies inside the circumcircle
        if (p[0] - self.center[0]) ** 2 + (p[1] - self.center[1]) ** 2 < self.radius ** 2:
            return True
        else:
            return False

    def plot(self):
        # plt.subplot(111).plot([self.__A[0], self.__B[0], self.__C[0], self.__A[0]], [self.__A[1], self.__B[1], self.__C[1], self.__A[1]])
        plt.subplot(111).add_artist(plt.Circle(self.center, self.radius, fill=False))
        # fig, ax = plt.subplots()
        # c = plt.Circle(self.center, self.radius, fill=False)
        # ax.add_artist(c)


