import numpy as np 
import pyvista as pv 
import matplotlib.pyplot as plt

data = pv.read("build/fm_0.vtk")
points = data.points
bounds = data.bounds
print(bounds)
pathlines = []
num_lines = 15 * 15 * 15
print(points.shape)
length = points.shape[0] / num_lines
print("length", length)

for n in range(num_lines):
    index = np.arange(n, num_lines * length, num_lines, dtype=np.int32)
    line = points[index, :]
    line = np.array(line)
    pathlines.append(line)

pathlines = np.array(pathlines)
print(pathlines.shape)
# print(pathlines[:, 0, :])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_xlim(bounds[0], bounds[1])
ax.set_ylim(bounds[2], bounds[3])
ax.set_zlim(bounds[4], bounds[5])

for i in range(num_lines):
    if i < 1000:
        ax.plot3D(pathlines[i, :, 0], pathlines[i, :, 1], pathlines[i, :, 2], c = 'blue')
        # ax.scatter3D(pathlines[i, 0, 0], pathlines[i, 0, 1], pathlines[i,0, 2], c='red')
plt.show()



    