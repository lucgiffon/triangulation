# code from http://codepad.org/KYyfy9CH
import math
import pylab as pl
import matplotlib.pyplot as plt
import numpy as np
import random
# Why math? Just 'cos


class CircumCircle:
    def __init__(self, A, B, C):
        self.__A = A
        self.__B = B
        self.__C = C

        self.center = self.circumcenter(self.__A, self.__B, self.__C)
        self.radius = self.circumradius(self.__A, self.__B, self.__C)

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

    def plot(self):
        # plt.subplot(111).plot([self.__A[0], self.__B[0], self.__C[0], self.__A[0]], [self.__A[1], self.__B[1], self.__C[1], self.__A[1]])
        plt.subplot(111).add_artist(plt.Circle(self.center, self.radius, fill=False))
        # fig, ax = plt.subplots()
        # c = plt.Circle(self.center, self.radius, fill=False)
        # ax.add_artist(c)


