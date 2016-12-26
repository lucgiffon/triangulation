import math
import random

import matplotlib.pyplot as plt
import numpy as np
# from matplotlib import collections as mc
import networkx as nx
import matplotlib.delaunay as triang

import circumcircle
from node import Node
from triangle import Triangle
from triangle_tracker import TriangleTracker


class Triangulator(nx.Graph):
    def __init__(self, x, y, bound=None):
        """
        :param bound: Indicate the maximum x and y of dots in the graph.
        """
        assert(len(x) == len(y))
        super().__init__()
        # the bound is the maximum absolute value on each axis where a node can be found
        self.__bound = bound
        # DAG which will keep track of triangles created and splitted
        self.__tracker = TriangleTracker()

        for i in range(len(x)):
            new_node = Node(x[i], y[i], i + 1)
            self.add_node(new_node)

        if self.__bound is None:
            self.__bound = max([max([math.fabs(n.x), math.fabs(n.y)]) for n in self.nodes_iter()])

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
            print("\t\t\t\tInvalid edge:", " ".join([str(n.x) + "," + str(n.y) for n in edge]))
            self.remove_edge(edge[0], edge[1])
            self.add_edge(new_node, point_k)
            self.__tracker.add_edge(illegality, Triangle(point_k, new_node, edge[0]))
            self.__tracker.add_edge(illegality, Triangle(point_k, new_node, edge[1]))
            self.legalize_edge(new_node, [point_k, edge[0]])
            self.legalize_edge(new_node, [point_k, edge[1]])

        else:
            print("\t\t\t\tValid edge: (%s, %s) - (%s, %s)" % (edge[0].x, edge[0].y, edge[1].x, edge[1].y))

    def is_edge_illegal(self, new_node, edge):

        # both the nodes of the candidate edge are negatives implies that the edge is one of the big triangle
        if edge[0].ident < 0 and edge[1].ident < 0:
            # we must keep the edges of the big triangle
            return None
        else:
            # get the triangles which contains the candidate edge
            triangles = self.__tracker.get_triangles_containing_edge(edge)
            for tri in triangles:
                if new_node in tri.vertices:
                    # if the new node is in the triangle, then this is not an interesting triangle: we are looking
                    # for the quad_point (point at the opposite side of the edge from the new_node).
                    continue
                else:
                    print("\t\t\tTested triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in tri]))
                # the quad point is the vertex of the triangle which is not at the boundaries of the candidate edge
                quad_point = set(tri.vertices).difference(set(edge)).pop()
                if edge[0].ident < 0 and edge[1].ident < 0:
                    print("\t\t\t\tValid edge: case I")
                    return None
                # if none of the 4 points is negative, then it is the normal case, we only check whether the quad_point
                # lies inside the circumcircle of the triangle tri
                elif min([edge[0].ident, edge[1].ident, new_node.ident, quad_point.ident]) > 0:
                    cc = circumcircle.CircumCircle((new_node.x, new_node.y), (edge[0].x, edge[0].y),
                                                   (edge[1].x, edge[1].y))
                    # plt.gca().set_aspect('equal', adjustable='box')
                    # cc.plot()
                    print("\t\t\tCircumcenter: %s - Radius: %s" % (cc.center, cc.radius))

                    # plt.plot([quad_point.x, edge[0].x, edge[1].x, new_node.x],
                    #          [quad_point.y, edge[0].y, edge[1].y, new_node.y], color='red')
                    # plt.plot([new_node.x, edge[0].x, edge[1].x, new_node.x],
                    #          [new_node.y, edge[0].y, edge[1].y, new_node.y], color='green')
                    # plt.scatter(new_node.x, new_node.y, marker="o")
                    # plt.scatter(quad_point.x, quad_point.y, marker="+")
                    # plt.scatter(tri[0].x, tri[0].y, marker="<")
                    # plt.scatter(tri[1].x, tri[1].y, marker="<")
                    # plt.scatter(tri[2].x, tri[2].y, marker="<")
                    # print(quad_point.x, "-", quad_point.y)
                    # print(cc.radius ** 2)
                    # print((quad_point.x - cc.center[0]) ** 2 + (quad_point.y - cc.center[1]) ** 2)
                    # plt.show()

                    # using the pythagore theorem, we test if the quad_point lies inside the circumcircle cc
                    if (quad_point.x - cc.center[0]) ** 2 + (quad_point.y - cc.center[1]) ** 2 < cc.radius ** 2:
                        # if so, the edge is illegal
                        print("\t\t\t\tInvalid edge: case II")
                        return tri
                    else:
                        print("\t\t\t\tValid edge: case II")
                        return None
                # if only one node of the edge is negative and nor the new_node nor the quad_point is,
                # then the edge is illegal: there shouldn't be an edge between a special node and a normal node
                elif (edge[0].ident < 0 or edge[1].ident < 0) and not (new_node.ident < 0 or quad_point.ident < 0):
                    print("\t\t\t\tInvalid edge: case III")
                    return tri

                # if no node of the edge is negative and either the new_node or the quad_point is negative
                # then the edge is illegal because we do not want edges between a node of the big triangle and a normal
                # node
                elif not (edge[0].ident < 0 or edge[1].ident < 0) and (new_node.ident < 0 or quad_point.ident < 0):
                    print("\t\t\t\tValid edge: case III")
                    return None

                # if there is one node of the big triangle and the quad_point negative, then, special case
                elif (edge[0].ident < 0 or edge[1].ident < 0) and (new_node.ident < 0 or quad_point.ident < 0):
                    # p1 and p2 are the two negative nodes
                    if edge[0].ident < 0:
                        p1 = edge[0]
                    else:
                        p1 = edge[1]
                    if new_node.ident < 0:
                        p2 = new_node
                    else:
                        p2 = quad_point

                    # if p1.ident < p2.ident:
                    #     print("\t\t\t\tValid edge: case IV")
                    #     return None
                    # else:
                    #     print("\t\t\t\tInvalid edge: case IV")
                    #     return tri
                    return None
                else:
                    print("new point: ", new_node.x, ", ", new_node.y)
                    print("edge: ", edge[0].x, ", ", edge[0].y, "  --  ", edge[1].x, ", ", edge[1].y)
                    print("quad point: ", quad_point.x, ", ", quad_point.y)

                    raise Exception("There is an unhandled case in the method Triangulator.is_edge_illegal()")
            return None

    def delaunay_triangulation(self):
        """
        Triangulate the nodes of the graph by adding edges between them.
        """

        # create the nodes corresponding to the big triangle
        p1 = Node(x=3 * self.__bound, y=0, ident=-1)
        p2 = Node(x=0, y=3 * self.__bound, ident=-2)
        p3 = Node(x=-3 * self.__bound, y=-3 * self.__bound, ident=-3)
        init_triangle = [p1, p2, p3]
        self.add_node(p1)
        self.add_node(p2)
        self.add_node(p3)
        self.add_edge(p1, p2)
        self.add_edge(p2, p3)
        self.add_edge(p1, p3)

        # # todo refaire les plot joliement avec un vrai subplot etc (animation)
        # plt.plot([p.x for p in init_triangle] + [p1.x], [p.y for p in init_triangle] + [p1.y])

        # the tracker is initialized with the big triangle
        self.__tracker.add_node(Triangle(*init_triangle))

        for new_node in set(self.nodes()).difference(set(init_triangle)):
            # TODO there will be multiple upper_triangle if the node is on an edge
            upper_triangle = self.__tracker.get_triangle_containing_node(new_node)
            print("New node: (%s, %s)" % (new_node.x, new_node.y))
            print("Upper triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in upper_triangle.vertices]))
            # input()
            for node_upper_triangle in upper_triangle.vertices:
                # at first, we create a new edge between the new_node and each vertex of the upper triangle
                # TODO destroy the edge under the new_node and create edge between the new_node and each vertex of
                #  both the upper triangles
                self.add_edge(node_upper_triangle, new_node)

            # trouver une autre option plus claire ou plus élégante, c'est nul, ça
            comb = [(0, 1), (1, 2), (0, 2)]

            # TODO check if there were two upper_triangle or only one
            for c in comb:
                # for each edge of the upper triangles (except the edge under the new_node), a new triangle is created
                new_tri = Triangle(new_node, upper_triangle[c[0]], upper_triangle[c[1]])
                print("\tNew triangle: ", " ".join([str(n.x) + "," + str(n.y) for n in new_tri.vertices]))
                # keep track of this new triangle
                self.__tracker.add_node(new_tri)
                # keep track of the ancestors of this new triangle
                self.__tracker.add_edge(upper_triangle, new_tri)
                # the new triangle may have induced illegal edges. Legalize them. The illegal edge should be the edge
                # at the opposite of the new_node on each triangle
                self.legalize_edge(new_node, [upper_triangle[c[0]], upper_triangle[c[1]]])

            self.plot()
            plt.show()


if __name__ == "__main__":

    N = 4
    BOUND = 100
    SHOW_EXAMPLE = False

    # plt.subplot(111).set_xlim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    # plt.subplot(111).set_ylim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    # plt.gca().set_aspect('equal', adjustable='box')

    # create triangulation
    # create dots to triangulate
    x, y = np.array(range(N)), np.array(range(N))
    x[0] = -50
    x[1] = 50
    x[2] = 50
    x[3] = -50

    y[0] = 50
    y[1] = 50
    y[2] = -50
    y[3] = -40
    # for i in range(N):
    #     x[i] = random.uniform(-BOUND, BOUND)
    #     y[i] = random.uniform(-BOUND, BOUND)

    triangulation = Triangulator(x, y, bound=BOUND)

    triangulation.delaunay_triangulation()
    triangulation.plot()

    plt.show()

    #############################################################
    # # # # # test code in order to compare the results # # # # #
    #############################################################
    if SHOW_EXAMPLE:
        cens,edg,tri,neig = triang.delaunay(x,y)
        for t in tri:
           # t[0], t[1], t[2] are the points indexes of the triangle
           t_i = [t[0], t[1], t[2], t[0]]
           # print(x[t_i])
           plt.plot(x[t_i], y[t_i])

        plt.plot(x, y, 'o')
        plt.show()
