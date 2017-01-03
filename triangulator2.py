"""
Delaunay triangulation.

Usage:
  triangulator NPOINTS BOUND [--show-example]

Arguments:
  NPOINTS   Number of points in the triangulation
  BOUND     Maximum absolute value of x and y value of points

Options:
  -h --help       Show this screen.
  --show-example  Show example of valid triangulation
"""

import math
import random

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import matplotlib.delaunay as triang
from docopt import docopt

import circumcircle
from node import Node
from triangle import Triangle
from triangle_tracker import TriangleTracker


class Triangulator:
    def __init__(self, bound):
        """
        Initialize the triangulator with x and y coordinates of all vertices of the graph.

        :param bound: Indicate the maximum absolute value of x and y of dots in the graph.
        """
        super().__init__()
        self.__bound = bound

        self.__nodes = [Node(*(self.__bound * 3 * np.array((-1, -1)))),
                        Node(*(self.__bound * 3 * np.array((1, -1)))),
                        Node(*(self.__bound * 3 * np.array((1, 1)))),
                        Node(*(self.__bound * 3 * np.array((-1, +1))))]

        # Create the first two triangles
        # The triangle is defined by its angles
        t1 = Triangle(self.__nodes[0], self.__nodes[1], self.__nodes[3])
        t2 = Triangle(self.__nodes[2], self.__nodes[3], self.__nodes[1])

        # Dict which keep track of triangles
        self.__tracker = dict()

        # Add the two triangles to the tracker
        # The tracker is indexed by the tuples representing the triangles
        # The values for each triangle are [opposite triangle toward edge 1, // edge 2, // edge 3]
        self.__tracker[t1] = [t2, None, None]
        self.__tracker[t2] = [t1, None, None]

        # https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
        # https://github.com/jmespadero/pyDelaunay2D/
        # https://en.wikipedia.org/wiki/Quadtree#Point_quadtree
        # https://github.com/karimbahgat/Pyqtree

    def plot(self):
        plt.plot([node.x for node in self.__nodes], [node.y for node in self.__nodes], 'o')
        for triangle in self.__tracker:
             plt.plot([triangle[0].x,
                       triangle[1].x,
                       triangle[2].x,
                       triangle[0].x],

                      [triangle[0].y,
                       triangle[1].y,
                       triangle[2].y,
                       triangle[0].y]
                      )

    def legalize_edge(self, new_node, edge):
        """
        Legalize the edge if it is illegal.

        :param new_node: The recently added node of the graph.
        :param edge: The edge that is maybe not legal
        """
        #!print("Processing edge: %s\t%s" % (edge[0], edge[1]))
        if self.is_edge_illegal(new_node, edge):
            # Get the triangles involved in the illegal edge
            triangles = self.__tracker.get_triangles_containing_edge(edge)

            assert len(triangles) == 2

            # Get the point_k. That is, the fourth point after the new node and the two nodes
            # of the edge
            opposite_vertexes = set(triangles[0].vertices + triangles[1].vertices).difference(set(edge))
            opposite_vertexes.remove(new_node)
            point_k = opposite_vertexes.pop()
            self.remove_edge(edge[0], edge[1])
            self.add_edge(new_node, point_k)

            new_triangle_1 = Triangle(point_k, new_node, edge[0])
            new_triangle_2 = Triangle(point_k, new_node, edge[1])

            for tri in triangles:
                self.__tracker.add_edge(tri, new_triangle_1)
                self.__tracker.add_edge(tri, new_triangle_2)

            self.legalize_edge(new_node, (point_k, edge[0]))
            self.legalize_edge(new_node, (point_k, edge[1]))
        else:
            pass

    def is_edge_illegal(self, new_node, edge):

        # get the triangles which contains the candidate edge
        triangles = self.__tracker.get_triangles_containing_edge(edge)
        assert len(triangles) == 1 or len(triangles) == 2

        for tri in triangles:
            if new_node in tri.vertices:
                # if the new node is in the triangle, then this is not an interesting
                # triangle: we are looking for the quad_point (point at the opposite side
                # of the edge from the new_node).
                continue
            else:
                pass
            quad_point = set(tri.vertices).difference(set(edge)).pop()
            if edge[0].ident < 0 and edge[1].ident < 0:
                return False
            # if none of the 4 points is negative, then it is the normal case, we only check
            # whether the quad_point
            # lies inside the circumcircle of the triangle tri
            # elif min([edge[0].ident, edge[1].ident, new_node.ident, quad_point.ident]) > 0:
            elif True:
                cc = circumcircle.CircumCircle((new_node.x, new_node.y), (edge[0].x, edge[0].y),
                                               (edge[1].x, edge[1].y))
                # using the pythagore theorem, we test if the quad_point lies inside the circumcircle cc
                if (quad_point.x - cc.center[0]) ** 2 + (quad_point.y - cc.center[1]) ** 2 < cc.radius ** 2:
                    # if so, the edge is illegal
                    return True
                else:
                    return False

            # if only one node of the edge is negative and nor the new_node nor the quad_point is,
            # then the edge is illegal: there shouldn't be an edge between a special node and a
            # normal node
            elif (edge[0].ident < 0 or edge[1].ident < 0) and not (new_node.ident < 0 or quad_point.ident < 0):
                return True
            # if no node of the edge is negative and either the new_node or the quad_point is negative
            # then the edge is legal because we do not want edges between a node of the big triangle
            # and a normal node
            elif not (edge[0].ident < 0 or edge[1].ident < 0) and (new_node.ident < 0 or quad_point.ident < 0):
                return False

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

                # lexicographic ordering
                if p1.ident < p2.ident:
                    return False
                else:
                    return True

            else:
                raise Exception("There is an unhandled case in the method Triangulator.is_edge_illegal()")
        # return False

    def add_node(self, node):
        node = np.asarray(node)
        node_index = len(self.__nodes)

        self.__nodes.append(node)

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

        # the tracker is initialized with the big triangle
        self.__tracker.add_node(Triangle(*init_triangle))

        for new_node in set(self.nodes()).difference(set(init_triangle)):
            upper_triangles = self.__tracker.get_triangle_containing_node(new_node)
            assert len(upper_triangles) == 1 or len(upper_triangles) == 2

            if len(upper_triangles) == 1:
                for node_upper_triangle in upper_triangles[0].vertices:
                    # create a new edge between the new_node and each vertex of the upper triangle
                    self.add_edge(node_upper_triangle, new_node)

                # Creating the new triangles
                new_tri_1 = Triangle(new_node,
                                     upper_triangles[0].vertices[0], upper_triangles[0].vertices[1])
                new_tri_2 = Triangle(new_node,
                                     upper_triangles[0].vertices[0], upper_triangles[0].vertices[2])
                new_tri_3 = Triangle(new_node,
                                     upper_triangles[0].vertices[1], upper_triangles[0].vertices[2])
                # keep track of these new triangles
                self.__tracker.add_node(new_tri_1)
                self.__tracker.add_node(new_tri_2)
                self.__tracker.add_node(new_tri_3)
                # keep track of the ancestors of this new triangle
                self.__tracker.add_edge(upper_triangles[0], new_tri_1)
                self.__tracker.add_edge(upper_triangles[0], new_tri_2)
                self.__tracker.add_edge(upper_triangles[0], new_tri_3)

                # the new triangle may have induced illegal edges. Legalize them.
                # The illegal edge should be the edge at the opposite of the new_node
                # on each triangle.
                self.legalize_edge(new_node,
                                   (upper_triangles[0].vertices[0], upper_triangles[0].vertices[1]))
                self.legalize_edge(new_node,
                                   (upper_triangles[0].vertices[0], upper_triangles[0].vertices[2]))
                self.legalize_edge(new_node,
                                   (upper_triangles[0].vertices[1], upper_triangles[0].vertices[2]))
            else:
                side_vertexes = tuple(set(upper_triangles[0].vertices).intersection(set(upper_triangles[1].vertices)))
                self.remove_edge(side_vertexes[0], side_vertexes[1])
                self.add_edge(new_node, side_vertexes[0])
                self.add_edge(new_node, side_vertexes[1])

                for triangle in upper_triangles:
                    opposite_vertex = set(triangle.vertices).difference(set(side_vertexes)).pop()
                    self.add_edge(opposite_vertex, new_node)

                    # Creating the new triangles
                    new_tri_1 = Triangle(new_node,
                                         side_vertexes[0], opposite_vertex)
                    new_tri_2 = Triangle(new_node,
                                         side_vertexes[1], opposite_vertex)
                    # keep track of these new triangles
                    self.__tracker.add_node(new_tri_1)
                    self.__tracker.add_node(new_tri_2)
                    # keep track of the ancestors of this new triangle
                    self.__tracker.add_edge(triangle, new_tri_1)
                    self.__tracker.add_edge(triangle, new_tri_2)

                    # the new triangle may have induced illegal edges. Legalize them.
                    # The illegal edge should be the edge at the opposite of the new_node
                    # on each triangle.
                    self.legalize_edge(new_node,
                                       (side_vertexes[0], opposite_vertex))
                    self.legalize_edge(new_node,
                                       (side_vertexes[1], opposite_vertex))


if __name__ == "__main__":

    arguments = docopt(__doc__)
    N = int(arguments["NPOINTS"])
    BOUND = float(arguments["BOUND"])
    SHOW_EXAMPLE = arguments["--show-example"]

    # create triangulation
    # create dots to triangulate
    x, y = np.array(range(N)), np.array(range(N))
    for i in range(N):
        x[i] = random.uniform(-BOUND, BOUND)
        y[i] = random.uniform(-BOUND, BOUND)

    triangulation = Triangulator(BOUND)

    plt.subplot(111).set_xlim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    plt.subplot(111).set_ylim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    plt.gca().set_aspect('equal', adjustable='box')

    triangulation.plot()
    plt.show()

    # triangulation.delaunay_triangulation()
    #

    # triangulation.plot()
    # plt.show()
    #
    # #############################################################
    # # # # # # test code in order to compare the results # # # # #
    # #############################################################
    # if SHOW_EXAMPLE:
    #     cens,edg,tri,neig = triang.delaunay(x,y)
    #     for t in tri:
    #        # t[0], t[1], t[2] are the points indexes of the triangle
    #        t_i = [t[0], t[1], t[2], t[0]]
    #        # #!print(x[t_i])
    #        plt.plot(x[t_i], y[t_i])
    #
    #     plt.subplot(111).set_xlim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    #     plt.subplot(111).set_ylim([-3 * BOUND * 1.5, 3 * BOUND * 1.5])
    #     plt.gca().set_aspect('equal', adjustable='box')
    #     plt.plot(x, y, 'o')
    #     plt.show()
