# import matplotlib.pyplot as plt
# import numpy as np
# N = 5
# x = np.random.rand(N)
# y = np.random.rand(N)
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

# import matplotlib.delaunay as triang
# import pylab
# import numpy
#
# # 10 random points (x,y) in the plane
# x,y =  numpy.array(numpy.random.standard_normal((2,10)))
# cens,edg,tri,neig = triang.delaunay(x,y)
#
# for t in tri:
#  # t[0], t[1], t[2] are the points indexes of the triangle
#  t_i = [t[0], t[1], t[2], t[0]]
#  print(t_i)
#  print(x[t_i], y[t_i])
#  input()
#  pylab.plot(x[t_i], y[t_i])
#  break
#
# pylab.plot(x, y, 'o')
# pylab.show()

# define set of point P (number N) -> graph (point = Node)
# define triangulation (= set of edges)
# define D, the structure containing triangles

class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

class Graph:
    pass