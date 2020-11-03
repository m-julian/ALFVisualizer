import pyvista as pv
import pandas as pd
import numpy as np
import meshio

df = pd.read_excel("methanol1000K_calculated_features.xlsx")

theta_data = np.array(df["θC1-H6"])
phi_data = np.array(df["φC1-H6"])

# Create a sphere
pi = np.pi
cos = np.cos
sin = np.sin

#DATA
x = sin(theta_data)*cos(phi_data)
y = sin(theta_data)*sin(phi_data)
z = cos(theta_data)
data_array = np.zeros((x.size,3))

for i in range(0,x.size):
    data_array[i][0] = x[i]
    data_array[i][1] = y[i]
    data_array[i][2] = z[i]

#making sphere and other things to plot
sphere = pv.Sphere(radius=1,center=(0,0,0))
# sphere_slices = sphere.slice_orthogonal(x=20,y=20,z=0)
sphere_slices = sphere.slice_orthogonal()
centre = np.zeros((1,3))
centre1 = pv.PolyData(centre)
z_axis = pv.Line(pointa=(0,0,0), pointb=(0,0,1))
z_axis_neg = pv.Line(pointa=(0,0,0), pointb=(0,0,-1))
x_axis = pv.Line(pointa=(0,0,0), pointb=(1,0,0))
y_axis = pv.Line(pointa=(0,0,0), pointb=(0,1,0))
point_cloud = pv.PolyData(data_array)

# xy_plane = np.zeros((x.size,2))
# for i in range(0,x.size):
#     xy_plane[i][0] = x[i]
#     xy_plane[i][1] = y[i]
# # xy_plane = pv.Plane(center=(0,0,0), direction=((xy_plane),0))
# xy_plane = pv.Plane(center=(0,0,0),direction=(0,0,1), i_size=2, j_size=1)

# xy_plane = pv.StructuredGrid(x, y, np.zeros((x.size,)))


# plotting
p = pv.Plotter(shape=(1,1))
p.subplot(0,0)
p.add_mesh(sphere,color='white', smooth_shading=True, opacity=0.50) #sphere
p.add_mesh(sphere_slices, color='green') #lines going around sphere
# p.add_mesh(xy_plane, color='green') #xy plane
p.add_mesh(centre1, color='red', point_size=10, render_points_as_spheres=True)
p.add_mesh(z_axis, color="blue", line_width=2) # z axis pisitive
p.add_mesh(z_axis_neg, color="blue", line_width=2) # z axis negative
p.add_mesh(x_axis, color="blue", line_width=2) # x axis pisitive
p.add_mesh(y_axis, color="blue", line_width=2) # y axis pisitive
p.add_mesh(point_cloud, color='orange', point_size=5, render_points_as_spheres=True) # theta and phi features
p.add_axes()


to_save = {'sphere': sphere, "sphere_slices":sphere_slices,
            'centre1':centre1, "z_axis":z_axis, "z_axis_neg":z_axis_neg, "x_axis": x_axis, "y_axis":y_axis, "point_cloud":point_cloud}
blocks = pv.MultiBlock(to_save)
merged = blocks.combine()

p.show(auto_close=False)

p.open_movie("test3.mp4")


x = np.arange(-10, 10, 0.25)
y = np.arange(-10, 10, 0.25)
x, y = np.meshgrid(x, y)
r = np.sqrt(x ** 2 + y ** 2)
z = np.sin(r)

# setup camera and close
p.show(auto_close=False)

pts = sphere.points.copy()

# Update Z and write a frame for each updated position
nframe = 15
for phase in np.linspace(0, 2 * np.pi, nframe + 1)[:nframe]:
    z = np.sin(r + phase)
    pts[:, -1] = z.ravel()
    p.update_coordinates(pts)
    p.update_scalars(z.ravel())
    p.write_frame()

# Close movie and delete object
p.close()