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

import random
import copy
import time
import itertools

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.delaunay as triang
from docopt import docopt
from pyqtree import Index

from node import Node
from triangle import Triangle


class Triangulator:
    def __init__(self, bound):
        """
        Initialize the triangulator with x and y coordinates of all vertices of the graph.

        :param bound: Indicate the maximum absolute value of x and y of dots in the graph.
        """
        super().__init__()
        self.__bound = bound

        p1 = Node(*(self.__bound * 3 * np.array((-1, -1))))
        p2 = Node(*(self.__bound * 3 * np.array((1, -1))))
        p3 = Node(*(self.__bound * 3 * np.array((1, 1))))
        p4 = Node(*(self.__bound * 3 * np.array((-1, +1))))
        self.__initializing_nodes = [p1, p2, p3, p4]

        self.__nodes = copy.copy(self.__initializing_nodes)

        # Create the first two ccw triangles
        # The triangle is defined by its angles
        t1 = Triangle(self.__nodes[0], self.__nodes[1], self.__nodes[3])
        t2 = Triangle(self.__nodes[2], self.__nodes[3], self.__nodes[1])

        # Dict which keep track of triangles
        self.__tracker = dict()
        self.__tracker_set = set()

        # Add the two triangles to the tracker
        # The tracker is indexed by the tuples representing the triangles
        # The values for each triangle are [opposite triangle toward edge 1, // edge 2, // edge 3]
        self.__tracker[t1] = [t2, None, None]
        self.__tracker[t2] = [t1, None, None]
        self.__tracker_set.add(t1)
        self.__tracker_set.add(t2)

        self.__spindex = Index(bbox=(*(self.__bound * 3 * np.array((-1, -1))),
                              *(self.__bound * 3 * np.array((1, 1)))))
        self.__spindex.insert(t1, t1.bbox)
        self.__spindex.insert(t2, t2.bbox)
        # https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
        # https://github.com/jmespadero/pyDelaunay2D/
        # https://en.wikipedia.org/wiki/Quadtree#Point_quadtree
        # https://github.com/karimbahgat/Pyqtree

        # plt.subplot(111).set_xlim([-BOUND * 3 * 1.1, BOUND * 3 * 1.1])
        # plt.subplot(111).set_ylim([-BOUND * 3 * 1.1, BOUND * 3 * 1.1])
        # ax = plt.gca()
        # ax.set_aspect('equal', adjustable='box')
        # ax.spines['right'].set_color('none')
        # ax.spines['top'].set_color('none')
        # ax.xaxis.set_ticks_position('bottom')
        # ax.spines['bottom'].set_position(('data', 0))
        # ax.yaxis.set_ticks_position('left')
        # ax.spines['left'].set_position(('data', 0))
        # plt.xticks([])
        # plt.yticks([])
        # self.plot_with_init_nodes()
        # plt.show()

    def plot(self, subplot=None):
        if subplot is None:
            plt.plot([node.x for node in self.__nodes if node not in self.__initializing_nodes],
                     [node.y for node in self.__nodes if node not in self.__initializing_nodes], 'o')
        else:
            plt.subplot(subplot).plot(
                [node.x for node in self.__nodes if node not in self.__initializing_nodes],
                [node.y for node in self.__nodes if node not in self.__initializing_nodes], 'o')
        for triangle in self.__tracker:
            if len(set(triangle.vertices).intersection(set(self.__initializing_nodes))) == 0:
                triangle.plot(subplot)

    def plot_with_init_nodes(self, subplot=None):
        if subplot is None:
            plt.plot([node.x for node in self.__nodes],
                     [node.y for node in self.__nodes], 'o')
        else:
            plt.subplot(subplot).plot([node.x for node in self.__nodes],
                                      [node.y for node in self.__nodes], 'o')
        for triangle in self.__tracker:
            triangle.plot(subplot)

    def add_node(self, node):
        node = Node(*node)

        self.__nodes.append(node)

        # search triangle whose circumcircle contains p
        candidate_bad_designed_triangles = set(self.__spindex.intersect((node.x, node.x, node.y, node.y)))
        # print(len(candidate_bad_designed_triangles.intersection(self.__tracker_set)), candidate_bad_designed_triangles.intersection(self.__tracker_set))
        # print(len(self.__tracker_set))
        bad_designed_triangles = []
        for triangle in candidate_bad_designed_triangles.intersection(self.__tracker_set):
        # for triangle in self.__tracker_set:
            if triangle.circumcircle_contains(node):
                bad_designed_triangles.append(triangle)

        # Find the external CCW boundary (star shape) of the bad triangles,
        # expressed as a list of edges (Node pairs) and the opposite
        # Triangle to each edge.

        # Initialization

        external_ccw_boundary = []
        bad_triangle = bad_designed_triangles[0]
        edge = 0
        while True:
            # get the opposite triangle around the current edge and for the current triangle
            opposite_triangle = self.__tracker[bad_triangle][edge]
            # check if the edge is on the boundary by looking if the opposite triangle is also a bad triangle
            if opposite_triangle not in bad_designed_triangles:
                # an entry of external_ccw_triangle is like: (edge, opposite_triangle)
                # edge being a pair of Nodes
                external_ccw_boundary.append(((bad_triangle[(edge+1) % 3],
                                               bad_triangle[(edge-1) % 3]),
                                              opposite_triangle))
                # look for next edge around the triangle
                edge = (edge+1) % 3

                # look if we went around the boundary of bad designed triangles
                # by comparing the first node of the first edge to the last node of the last edge
                if external_ccw_boundary[0][0][0] == external_ccw_boundary[-1][0][1]:
                    break
            else:
                # Move to next CCW edge in opposite triangle
                edge = (self.__tracker[opposite_triangle].index(bad_triangle) + 1) % 3
                bad_triangle = opposite_triangle

        # Remove the bad_designed_triangles
        for triangle in bad_designed_triangles:
            del self.__tracker[triangle]
            self.__tracker_set.remove(triangle)

        # Retriangulate the hole left by bad_triangles
        new_triangles = []
        for (edge, opposite_triangle) in external_ccw_boundary:
            new_triangle = Triangle(node,
                             edge[0],
                             edge[1])
            self.__tracker[new_triangle] = [opposite_triangle, None, None]
            self.__tracker_set.add(new_triangle)
            self.__spindex.insert(new_triangle, new_triangle.bbox)

            # If the new_triangle is not at the boundary of the triangulation
            if opposite_triangle is not None:
                # set the new_triangle as neighbour of the opposite_triangle
                for i, neighbour in enumerate(self.__tracker[opposite_triangle]):
                    if neighbour is not None and neighbour.contains_edge(edge):
                        # change link to use our new_triangle
                        self.__tracker[opposite_triangle][i] = new_triangle

            new_triangles.append(new_triangle)

        N = len(new_triangles)
        for i, triangle in enumerate(new_triangles):
            self.__tracker[triangle][1] = new_triangles[(i+1) % N]
            self.__tracker[triangle][2] = new_triangles[(i-1) % N]

        # self.plot()
        # plt.show()


def generate_list_of_points(size, min_, max_):
    points = {(random.randint(min_, max_), random.randint(min_, max_)) for i in range(N)}
    while len(points) < size:
        points |= {(random.randint(min_, max_), random.randint(min_, max_))}
    return points
    # points = list(list(x) for x in points)

if __name__ == "__main__":

    arguments = docopt(__doc__)
    N = int(arguments["NPOINTS"])
    BOUND = float(arguments["BOUND"])
    SHOW_EXAMPLE = arguments["--show-example"]

    if int(BOUND) ** 2 <= N:
        exit("Too mutch point asked in the specified area. Please, select a number of points lesser than the square of your bound.")

    # create triangulation
    triangulation = Triangulator(BOUND)

    # create dots to triangulate
    X, Y = zip(*generate_list_of_points(N, -BOUND, BOUND))

    start = time.time()
    for i in range(N):
        triangulation.add_node((X[i], Y[i]))
    stop = time.time()
    print("execution_time", stop - start)

    plt.subplot(111).set_xlim([-BOUND * 1.1, BOUND * 1.1])
    plt.subplot(111).set_ylim([-BOUND * 1.1, BOUND * 1.1])
    plt.gca().set_aspect('equal', adjustable='box')
    triangulation.plot()
    plt.show()

    #############################################################
    # # # # # test code in order to compare the results # # # # #
    #############################################################
    if SHOW_EXAMPLE:
        cens,edg,tri,neig = triang.delaunay(X, Y)
        for t in tri:
           # t[0], t[1], t[2] are the points indexes of the triangle
           t_i = [t[0], t[1], t[2], t[0]]
           # #!print(X[t_i])
           plt.plot(X[t_i], Y[t_i])

        plt.subplot(111).set_xlim([-BOUND * 1.5, BOUND * 1.5])
        plt.subplot(111).set_ylim([-BOUND * 1.5, BOUND * 1.5])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(X, Y, 'o')
        plt.show()
