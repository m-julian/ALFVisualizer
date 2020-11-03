import pyvista as pv
import pandas as pd
import numpy as np
import meshio
import calculate_alf

from glob import glob

xyz_files = sorted(glob("*.xyz"))
if len(xyz_files) == 1:
    xyz_file = xyz_files[0]
else:
    print("Select xyz file to evaluate:")
    print ("")
    for i in range(len(xyz_files)):
        print (f"[{i+1}] --> ", xyz_files[i])
    xyz_file = xyz_files[int(input())-1]

all_atom_features, atom_names = calculate_alf.features_and_atom_names(xyz_file)
# all_atom_features are [atom][point][feature]

print(all_atom_features[0,:,:])

center = np.array([0,0,0])
center_point = pv.PolyData(center)

# C1 bond 1
bond1 = all_atom_features[5,:,[0]].T
print(bond1)
# min_bond1 = min(bond1[:,0])
# max_bond1 = max(bond1[:,0])
# bond1 = np.divide(np.subtract(bond1, min_bond1),(max_bond1 - min_bond1))
z = np.zeros((bond1.shape[0], 2), dtype=bond1.dtype)
xy_bond1 = np.concatenate((bond1,z), axis=1)
bond1_cloud = pv.PolyData(xy_bond1)

# C1 bond 2
bond2 = all_atom_features[0,:,[1]].T
angle12 = all_atom_features[0,:,[2]].T
x_bond2 = np.multiply(np.cos(angle12), bond2)
y_bond2 = np.multiply(np.sin(angle12), bond2)
z_bond2 = np.zeros((bond2.shape[0], 1), dtype=bond2.dtype)
xy_bond2 = np.concatenate((x_bond2,y_bond2,z_bond2), axis=1)
bond2_cloud = pv.PolyData(xy_bond2)

r_data = all_atom_features[0,:,[3]].T
theta_data = all_atom_features[0,:,[4]].T
phi_data = all_atom_features[0,:,[5]].T
xx = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.cos(phi_data))
yy = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.sin(phi_data))
zz = np.multiply(r_data, np.cos(theta_data))
print(xx)

r3 = np.concatenate((xx,yy,zz), axis=1)
r3_cloud = pv.PolyData(r3)

r_data = all_atom_features[0,:,[6]].T
theta_data = all_atom_features[0,:,[7]].T
phi_data = all_atom_features[0,:,[8]].T
xx = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.cos(phi_data))
yy = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.sin(phi_data))
zz = np.multiply(r_data, np.cos(theta_data))
print(xx)

r4 = np.concatenate((xx,yy,zz), axis=1)
r4_cloud = pv.PolyData(r4)

r_data = all_atom_features[0,:,[9]].T
theta_data = all_atom_features[0,:,[10]].T
phi_data = all_atom_features[0,:,[11]].T
xx = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.cos(phi_data))
yy = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.sin(phi_data))
zz = np.multiply(r_data, np.cos(theta_data))
print(xx)

r5 = np.concatenate((xx,yy,zz), axis=1)
r5_cloud = pv.PolyData(r5)


# for r_col, theta_col, phi_col in zip(r_cols,theta_cols,phi_cols):

#     random_color = tuple([random.random() for i in range(0,3)])

#     r_data = np.array(df[r_col])
#     theta_data = np.array(df[theta_col])
#     phi_data = np.array(df[phi_col])

#     xx = r_data*sin(theta_data)*cos(phi_data)
#     yy = r_data*sin(theta_data)*sin(phi_data)
#     zz = r_data*cos(theta_data)

#     plotter.points3d(xx, yy, zz, scale_factor=0.03, color=(random_color))


plotter = pv.Plotter()
plotter.add_mesh(center_point, color='red', point_size=20, render_points_as_spheres=True)
plotter.add_mesh(bond1_cloud, color='maroon', point_size=10.,
                 render_points_as_spheres=True)
plotter.add_mesh(bond2_cloud, color='blue', point_size=10.,
                 render_points_as_spheres=True)
plotter.add_mesh(r3_cloud, color='#b14730', point_size=10.,
                 render_points_as_spheres=True)
plotter.add_mesh(r4_cloud, color='#FFFFFF', point_size=10.,
                 render_points_as_spheres=True)
plotter.add_mesh(r5_cloud, color='grey', point_size=10.,
                 render_points_as_spheres=True)
plotter.show_grid()
plotter.show()


# # x = np.array([[0,0,0],[1,0,0], [1,1,0]])
# # print(x.shape)
# # print(x)
# # point_cloud = pv.PolyData(x)

# # # Create random XYZ points
# points = np.random.rand(100, 3)
# # print(points.shape)
# # # Make PolyData
# point_cloud = pv.PolyData(points)
# # point_cloud.plot(eye_dome_lighting=True)

# points2 = -np.random.rand(100, 3)
# point_cloud2 = pv.PolyData(points2)

# plotter = pv.Plotter()
# plotter.add_mesh(center_point, color='red', point_size=20, render_points_as_spheres=True)
# plotter.add_mesh(point_cloud, color='maroon', point_size=10.,
#                  render_points_as_spheres=True)
# plotter.add_mesh(point_cloud2, color='maroon', point_size=10.,
#                  render_points_as_spheres=True)
# plotter.show_grid()
# plotter.show()