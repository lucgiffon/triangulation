from math import sqrt
from node import Node
from circumcircle import CircumCircle
import matplotlib.pyplot as plt

class Triangle:
    def __init__(self, p1, p2, p3):
        self.vertices = [p1, p2, p3]

        self.__circumcircle = CircumCircle((p1.x, p1.y),
                                           (p2.x, p2.y),
                                           (p3.x, p3.y))

        self.bbox = (self.__circumcircle.center[0] - self.__circumcircle.radius,
                     self.__circumcircle.center[1] - self.__circumcircle.radius,
                     self.__circumcircle.center[0] + self.__circumcircle.radius,
                     self.__circumcircle.center[1] + self.__circumcircle.radius)

    def __getitem__(self, item):
        return self.vertices[item]

    def __iter__(self):
        return iter(self.vertices)

    def circumcircle_contains(self, p):
        return self.__circumcircle.contains((p.x, p.y))

    def contains(self, p):

        p1 = self.vertices[0]
        p2 = self.vertices[1]
        p3 = self.vertices[2]

        alpha = ((p2.y - p3.y) * (p.x - p3.x) + (p3.x - p2.x) * (p.y - p3.y)) / ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y))

        beta = ((p3.y - p1.y) * (p.x - p3.x) + (p1.x - p3.x) * (p.y - p3.y)) / ((p2.y - p3.y) * (p1.x - p3.x) + (p3.x - p2.x) * (p1.y - p3.y))

        gamma = 1.0 - alpha - beta

        return alpha >= 0 and beta >= 0 and gamma >= 0

    def contains_edge(self, edge):
        return set(edge).issubset(set(self.vertices))

    def plot(self):
        plt.plot([v.x for v in self.vertices + [self.vertices[0]]],
                 [v.y for v in self.vertices + [self.vertices[0]]])

    def __repr__(self):
        return "\t".join(["(%s, %s)" % (v.x, v.y) for v in self.vertices])