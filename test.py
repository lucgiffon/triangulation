import math
import random

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import collections as mc
import networkx as nx
import numpy

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
        plt.plot([node.x for node in self.nodes()], [node.y for node in self.nodes()], 'o')
        for e in self.edges():
            plt.plot([e[0].x, e[1].x], [e[0].y, e[1].y])

    def test(self):
        for node in self.nodes_iter():
            for node2 in self.nodes_iter():
                self.add_edge(node, node2)

    def legalize_edge(self, new_node, edge):
        print("\t\tEdge: ", " ".join([str(n.x) + "," + str(n.y) for n in edge]))
        illegality = self.is_edge_illegal(new_node, edge)
        if illegality is not None:
            # replace edge by new_node--k
            point_k = set(illegality.vertices).difference(set(edge)).pop()
            self.remove_edge(edge[0], edge[1])
            self.add_edge(new_node, point_k)
            self.__tracker.add_edge(illegality, Triangle(point_k, new_node, edge[0]))
            self.__tracker.add_edge(illegality, Triangle(point_k, new_node, edge[1]))
            self.legalize_edge(new_node, [point_k, edge[0]])
            self.legalize_edge(new_node, [point_k, edge[1]])

    def is_edge_illegal(self, new_node, edge):
        triangles = self.__tracker.get_triangles_containing_edge(edge)
        for tri in triangles:
            print("\t\tTested triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in tri]))
            if new_node in tri.vertices:
                pass
            else:

                cc = circumcircle.CircumCircle((new_node.x, new_node.y), (edge[0].x, edge[0].y), (edge[1].x, edge[1].y))
                plt.gca().set_aspect('equal', adjustable='box')
                # cc.plot()
                print("\t\t\tCircumcenter: %s - Radius: %s" % (cc.center, cc.radius))
                quad_point = set(tri.vertices).difference(set(edge)).pop()
                # plt.plot([quad_point.x, edge[0].x, edge[1].x, new_node.x],
                #          [quad_point.y, edge[0].y, edge[1].y, new_node.y], color='red')
                # plt.plot([new_node.x, edge[0].x, edge[1].x, new_node.x],
                #          [new_node.y, edge[0].y, edge[1].y, new_node.y], color='green')
                # plt.scatter(new_node.x, new_node.y, marker="o")
                # plt.scatter(quad_point.x, quad_point.y, marker="+")
                # plt.scatter(tri[0].x, tri[0].y, marker="<")
                # plt.scatter(tri[1].x, tri[1].y, marker="<")
                # plt.scatter(tri[2].x, tri[2].y, marker="<")
                print(quad_point.x, "-", quad_point.y)
                print(cc.radius ** 2)
                print((quad_point.x - cc.center[0]) ** 2 + (quad_point.y - cc.center[1]) ** 2)
                # plt.show()
                if (quad_point.x - cc.center[0]) ** 2 + (quad_point.y - cc.center[1]) ** 2 < cc.radius ** 2:
                    return tri

        return None

    def delaunay_triangulation(self):

        if self.__bound is None:
            self.__bound = max([max([math.fabs(n.x), math.fabs(n.y)]) for n in self.nodes_iter()])

        p1 = Node(x=3 * self.__bound, y=0)
        p2 = Node(x=0, y=3 * self.__bound)
        p3 = Node(x=-3 * self.__bound, y=-3 * self.__bound)
        init_triangle = [p1, p2, p3]

        plt.plot([p.x for p in init_triangle] + [p1.x], [p.y for p in init_triangle] + [p1.y])

        self.__tracker.add_node(Triangle(*init_triangle))

        for new_node in self.nodes():
            upper_triangle = self.__tracker.get_triangle_containing_node(new_node)
            print("Upper triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in upper_triangle]))
            for node_upper_triangle in upper_triangle.vertices:
                self.add_edge(node_upper_triangle, new_node)

            comb = [(0, 1), (1, 2), (0, 2)]

            for c in comb:
                new_tri = Triangle(new_node, upper_triangle[c[0]], upper_triangle[c[1]])
                print("\tNew triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in new_tri.vertices]))
                self.__tracker.add_node(new_tri)
                self.__tracker.add_edge(upper_triangle, new_tri)
                self.legalize_edge(new_node, [upper_triangle[c[0]], upper_triangle[c[1]]])


class Triangle:
    def __init__(self, p1, p2, p3):
        self.vertices = [p1, p2, p3]

    def __getitem__(self, item):
        return self.vertices[item]

    def __iter__(self):
        return iter(self.vertices)

    def contains(self, node):
        def sign(p1, p2, p3):
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)

        # verifier que ca marche pour les noeuds qui sont sur une arrete
        b1 = sign(node, self.vertices[0], self.vertices[1]) < 0.0
        b2 = sign(node, self.vertices[1], self.vertices[2]) < 0.0
        b3 = sign(node, self.vertices[2], self.vertices[0]) < 0.0

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
            if set(edge).issubset(set(triangle.vertices)):
                triangles.append(triangle)
        return triangles


# angle = random.randint(1, 360) * math.pi * 2
# for i in range(0, 100):
#     pl.plot([math.cos(angle * i) * 50, math.sin(angle * i) * 50], 'o')


N = 10
BOUND = 100

plt.subplot(111).set_xlim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
plt.subplot(111).set_ylim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
plt.gca().set_aspect('equal', adjustable='box')

triangulation = Graph(bound=BOUND)
# x,y = numpy.array(numpy.random.standard_normal((2,N)))
x, y = np.array(range(N)), np.array(range(N))

for i in range(N):
    x[i] = random.uniform(-BOUND, BOUND)
    y[i] = random.uniform(-BOUND, BOUND)
    new_node = Node(x[i], y[i])
    triangulation.add_node(new_node)

triangulation.delaunay_triangulation()
triangulation.plot()

plt.show()

import matplotlib.delaunay as triang
import pylab

# 10 random points (x,y) in the plane

cens,edg,tri,neig = triang.delaunay(x,y)
for t in tri:
   # t[0], t[1], t[2] are the points indexes of the triangle
   t_i = [t[0], t[1], t[2], t[0]]
   # print(x[t_i])
   plt.plot(x[t_i], y[t_i])

plt.plot(x, y, 'o')
plt.show()

# define set of point P (number N) -> graph (point = Node)
# define triangulation (= set of edges)
# define D, the structure containing triangles