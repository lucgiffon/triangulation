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

        b1 = sign(node, self.vertices[0], self.vertices[1]) <= 0.0
        b2 = sign(node, self.vertices[1], self.vertices[2]) <= 0.0
        b3 = sign(node, self.vertices[2], self.vertices[0]) <= 0.0

        return (b1 == b2) and (b2 == b3)

    def __repr__(self):
        return "\t".join(["(%s, %s)" % (v.x, v.y) for v in self.vertices])