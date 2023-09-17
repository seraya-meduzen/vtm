import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np
from skimage import measure

x = np.linspace(-3,3, 400)
y = np.linspace(-3,3, 200)
z = np.linspace(-3,3, 200)

X, Y, Z =  np.meshgrid(x, y, z)

def f(x,y,z):
    return (x ** 2 + (9 / 4) * y ** 2 + z ** 2 - 1) ** 3  - x ** 2 * z ** 3 - (9 / 200) * y ** 2 * z ** 3

verts, faces, normals, values = measure.marching_cubes(f(X,Y,Z), 0)

fig = plt.figure(figsize=(30, 20))
ax = fig.add_subplot(111, projection='3d')

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["navy", "maroon", "darkviolet"])
ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2], cmap=cmap, lw=1)

ax.view_init(5, 0)

ax.grid(False)
ax.set_axis_off()

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

enc = b'\xd0\x9d\xd0\xb0\xd1\x81\xd1\x82\xd1\x8c), \xd0\xb3\xd0\xbe \xd0\xb2\xd1\x81\xd1\x82\xd1\x80\xd0\xb5\xd1\x82\xd0\xb8\xd0\xbc\xd1\x81\xd1\x8f))00)'

ax.set_title(enc.decode('utf8', 'strict'), fontsize=25)
plt.show()
