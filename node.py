class Node:
    def __init__(self, x, y, ident=None):
        """
        Node of the graph constitued by x and y coordinates and an index:

            - If the index is positive, the Node is legit
            - Else, this is a special node

        :param x: x coordinate of node
        :param y: y coordinate of node
        :param ident: index of the node
        """
        self.x = x
        self.y = y
        self.ident = ident

    # shouldn't be usefull
    # todo remove this
    # def __eq__(self, other):
    #     return (self.x == other.x) and (self.y == other.y)

    def __repr__(self):
        return "(%s, %s)" % (self.x, self.y)