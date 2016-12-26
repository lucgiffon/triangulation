import networkx as nx

class TriangleTracker(nx.DiGraph):
    """
    DAG tracking the appearance and splitting of triangles while triangulating.

    The leaves correspond to the triangles of the current triangulation.
    Nodes inside the dag refer to the triangle that have been created then splitted into new triangles.
    Edges refer to the splitting of a Triangle into new triangles.
    """

    def get_triangle_containing_node(self, node, p_node=None):
        """
        Given a node (a vertex in the graph), return the triangle of the current triangulation which contain that node.

        The analysis is done top down instead of brwosing the leaves directly because there is at max log(n) steps in the
        tree but 2n - 2 -k leaves (n being the number of nodes and k the number of edges on the convex hull).

        :param node: The node being looked for
        :param p_node: Used for recursion
        :return: The triangle containing the vertex
        """
        # TODO this should return a list: if the node is on an edge, there can be 2 triangles containing the node.
        if p_node is None:
            # get the triangle at the top of the tree: the big triangle.
            triangles = (x for x in self.nodes_iter() if self.in_degree(x) == 0)
        else:
            # get the next triangles for the topdown analysis.
            triangles = self.successors_iter(p_node)

        for triangle in triangles:
            if triangle.contains(node):
                # We are looking for the tighest triangles, if the triangle has successors,
                #  it means that it have been splitted
                if self.successors(triangle):
                    # recursion
                    return self.get_triangle_containing_node(node, p_node=triangle)
                else:
                    # there is no successor: the current triangle is the tightest triangle found
                    return triangle

    def get_triangles_containing_edge(self, edge):
        """
        Return the triangles containing the edge.

        :param edge: The edge being looked for
        :return: List of the triangles containing the edge
        """
        # found triangle are stored in list because there is one or two triangles containing a given edge
        triangles = []
        for triangle in (x for x in self.nodes_iter() if self.out_degree(x) == 0):
            if set(edge).issubset(set(triangle.vertices)):
                triangles.append(triangle)
        return triangles