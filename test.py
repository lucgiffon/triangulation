import math
import random

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
from matplotlib import collections as mc
import networkx as nx

import circumcircle

class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Graph(nx.Graph):
    def __init__(self, bound=None):
        super().__init__()
        self.__bound = bound
        self.__tracker = TriangleTracker()


    def plot(self):
        pl.plot([node.x for node in self.nodes()], [node.y for node in self.nodes()], 'o')
        for e in self.edges():
            pl.plot([e[0].x, e[1].x], [e[0].y, e[1].y])
        pl.show()

    def test(self):
        for node in self.nodes_iter():
            for node2 in self.nodes_iter():
                self.add_edge(node, node2)

    def legalize_edge(self, new_node, edge):
        print("\t\tEdge: ", " ".join([str(n.x) + "," + str(n.y) for n in edge]))
        if self.is_edge_illegal(new_node, edge):
            # replace edge by new_node--k
            #self.legalize_edge(new_node, edge[0]--k)
            #self.legalize_edge(new_node, edge[1]--k)
            pass

    def is_edge_illegal(self, new_node, edge):
        triangles = self.__tracker.get_triangles_containing_edge(edge)
        for tri in triangles:
            print("\t\tTested triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in tri]))
            circumcenter = circumcircle.circumcenter(*[(n.x, n.y) for n in tri])
            radius = circumcircle.radius(*[(n.x, n.y) for n in tri])
            print(circumcenter)
            angle = random.randint(1, 360) * math.pi * 2
            for i in range(0, 100):
                pl.plot([math.cos(angle * i) * radius, math.sin(angle * i) * radius])

        pass

    def delaunay_triangulation(self):

        if self.__bound is None:
            self.__bound = max([max([math.fabs(n.x), math.fabs(n.y)]) for n in self.nodes_iter()])

        p1 = Node(x=3 * self.__bound, y=0)
        p2 = Node(x=0, y=3 * self.__bound)
        p3 = Node(x=-3 * self.__bound, y=-3 * self.__bound)
        init_triangle = [p1, p2, p3]

        self.__tracker.add_node(Triangle(*init_triangle))

        for node in self.nodes():
            upper_triangle = self.__tracker.get_triangle_containing_node(node)
            print("Upper triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in upper_triangle]))
            for node_triangle in upper_triangle:
                self.add_edge(node_triangle, node)

            comb = [(0, 1), (1, 2), (0, 2)]

            for c in comb:
                tri = Triangle(node, upper_triangle[c[0]], upper_triangle[c[1]])
                print("\tNew triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in tri]))
                self.__tracker.add_node(tri)
                self.__tracker.add_edge(upper_triangle, tri)
                self.legalize_edge(node, [upper_triangle[c[0]], upper_triangle[c[1]]])


class Triangle:
    def __init__(self, p1, p2, p3):
        self.nodes = [p1, p2, p3]

    def __getitem__(self, item):
        return self.nodes[item]

    def __iter(self):
        return self.nodes

    def contains(self, node):
        def sign(p1, p2, p3):
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)

        # verifier que ca marche pour les noeuds qui sont sur une arrete
        b1 = sign(node, self.nodes[0], self.nodes[1]) < 0.0
        b2 = sign(node, self.nodes[1], self.nodes[2]) < 0.0
        b3 = sign(node, self.nodes[2], self.nodes[0]) < 0.0

        return (b1 == b2) and (b2 == b3)


class TriangleTracker(nx.DiGraph):

    def get_triangle_containing_node(self, node, p_node=None):
        if p_node is None:
            triangles = (x for x in self.nodes_iter() if self.in_degree(x)==0)
        else:
            triangles = self.successors_iter(p_node)

        for triangle in triangles:
            if triangle.contains(node):
                if self.successors(triangle):
                    return self.get_triangle_containing_node(node, p_node=triangle)
                else:
                    return triangle

    def get_triangles_containing_edge(self, edge):
        triangles = []
        for triangle in (x for x in self.nodes_iter() if self.out_degree(x)==0):
            if set(edge).issubset(set(triangle)):
                triangles.append(triangle)
        return triangles



N = 5
BOUND = 100
triangulation = Graph(bound=BOUND)
for i in range(N):
    new_node = Node(random.uniform(-BOUND, BOUND), random.uniform(-BOUND, BOUND))
    triangulation.add_node(new_node)

triangulation.delaunay_triangulation()
triangulation.plot()



# # plt.scatter(x, y)
# # plt.show()
# import pylab as pl
# from matplotlib import collections  as mc
# lines = []
# for i in range(N):
#     for j in range(N):
#             lines.append([(x[i], y[i]), (x[j], y[j])])
# lc = mc.LineCollection(lines, linewidths=2)
# fig, ax = pl.subplots()
# ax.add_collection(lc)
# ax.autoscale()
# ax.margins(0.1)
# plt.show()

import matplotlib.delaunay as triang
import pylab
import numpy

# 10 random points (x,y) in the plane
x,y =  numpy.array(numpy.random.standard_normal((2,10)))
cens,edg,tri,neig = triang.delaunay(x,y)
for t in tri:
   # t[0], t[1], t[2] are the points indexes of the triangle
   t_i = [t[0], t[1], t[2], t[0]]
   # print(x[t_i])
   pylab.plot(x[t_i], y[t_i])

pylab.plot(x, y, 'o')
pylab.show()

# define set of point P (number N) -> graph (point = Node)
# define triangulation (= set of edges)
# define D, the structure containing triangles