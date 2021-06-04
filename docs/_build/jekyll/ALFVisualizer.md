---
date: '2021-02-20T23:34:43.453Z'
docname: ALFVisualizer
images: {}
path: /alf-visualizer
title: ALFVisualizer module
---

# ALFVisualizer module


### class ALFVisualizer.Atom(coordinate_line)
Bases: `object`

Deals with 1 Atom / 1 Coordinate Line, self._atoms containts items of class Atom,
they are printed nicely because of the __str__ method of Atom class


#### C_1k()

#### C_2k()

#### C_3k(C_1k, C_2k)

#### add_alf_atom(atom)

#### property alf()

#### property alf_nums()

#### ang2bohr( = 1.88971616463)

#### angle(atom1, atom2)

#### property atom()

#### property atom_coordinates()

#### property atom_num()

#### property atom_num_coordinates()

#### property bonds()

#### calculate_features(atoms, unit='bohr')

#### property coordinates()

#### property coordinates_string()

#### counter( = count(1))

#### dist(other)

#### get_priorty(level)

#### property mass()

#### property name()

#### property num()

#### property priority()

#### property radius()

#### read_input(coordinate_line)

#### reset_level()

#### set_bond(jatom)

#### to_angstroms()

#### to_bohr()

#### property type()

#### xdiff(other)

#### ydiff(other)

#### zdiff(other)

### class ALFVisualizer.Atoms(atoms=None)
Bases: `object`

Deals with 1 trajectory timestep / all atoms in a timestep


#### ALF( = [])

#### PRIORITIES( = [])

#### add(atom)
adds an instance of Atom to _atoms


#### property alf()

#### property atoms()

#### calculate_alf()

#### calculate_features()
calculates features for atoms in one timestep


#### connect(iatom, jatom)

#### property connectivity()

#### property empty()

#### property features()

#### property masses()

#### property max_priority()

#### property nfeatures()

#### property priority()

#### set_alf()

#### to_angstroms()

#### to_bohr()

### class ALFVisualizer.OpenXYZFile()
Bases: `PyQt5.QtWidgets.QWidget`

Open File Dialog to select XYZ file


#### open_file_dialog()

#### open_pyvista_box(xyz_file)

### class ALFVisualizer.Trajectory(fname, read=False)
Bases: `object`

Reads trajectory file


#### add(atoms)
each item in self._trajectory is of Atoms class,
printing self._trajectory prints out based on Atoms __str__ method


#### property atom_names()

#### property features()

#### property nfeatures()

#### property priorities()

#### read()
read in trajectory file (typically .xyz)


### class ALFVisualizer.VisualizationWindow(all_atom_dict, atom_names)
Bases: `PyQt5.QtWidgets.QMainWindow`

handles GUI and connects user commands with what to plot on pyvista plot
see [https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog](https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog) videos for info on using Qt with python


#### change_atom_color()
opens up color dialog and lets user select a new color for a particular atom. GUI is automatically frozen until user closes the color dialog box.


#### default_or_random_atom_colors()
updates atom colors depending on default_atom_colors_checkbox. This checkbox is initialized to unchecked state, so
random colors are used (this is initial state). If checkbox is checked, then default colors for atoms,
i.e. oxygen:red , hydrogen:white, etc. are used


#### grid_status()
show or remove grid on pyvista plot depending on grid checkbox, updates atom data to plot


#### plot_updated_data()
plots all the data after all the checkboxes/sliders/colors etc. have been processed


#### points_slider_reset()
reset slider when changing to a different central atom


#### slider_update_plotted_data()
Makes slices of the data before making the data a pyvista MultiBlock, see update_noncentral_atoms_data method


#### untick_all_noncentral_atoms()
method called after clicking the Remove All Plotted Atoms button - this unticks all noncentral atoms


#### update_atom_color_box_buttons()
updates atoms that are in the color box, depending on which central atom is chosen and also which
checkboxes are ticked in the checkbox widget.


#### update_atom_data_and_plot()

#### update_central_atom_and_plot()

#### update_central_atom_data()
method used to update the central ALF atom, depending on selected atom in combo box


#### update_checkboxes_widget()
Used to dynamically generate the non-central atom checkboxes. They can be used to plot individual noncentral atoms instead of all noncentral atoms.


#### update_checked_atoms()
Method that keeps track of which checkboxes for noncentral atoms are checked. ONLY atoms that are checked are plotted


#### update_noncentral_atoms_data()
updates non-central atoms data to be plotted. This method gives all the non-central atoms and manipulates the data to be plotted if the slider has been
changed from its inital value. It also makes a pyvista MultiBlock object which is then used for plotting the non-central atoms and their colors.


### class ALFVisualizer.VisualizationWindowDecorators()
Bases: `object`

Decorators for UI


#### static clear_plot_use_grid(original_method)
remove or show grid for pyvista


### class ALFVisualizer.XYZArrays(all_atom_features)
Bases: `object`

Class for converting to Cartesian space. 
Creates a 3D array for each atom on which the ALF is centered (1 4D array total). Each 2D array in the 3D array
consists of N_pointsx3 (because each point has x,y,z coords in 3D space) matrices, where each matrix contains
the xyz coordinates of every atom that is NOT the atom on which the ALF is centered.


#### get_polar_atom_3d_array(one_atom)
Input: Takes in one atom feature matrix. 
Every three features (after the first three features that define the xy plane) have a radius, theta, and
phi component that defines where an atom is in xyz coordinates
This method will return a 3D array of shape n_remaining_atoms, n_points, 3, wher 3 is because x,y,z coordinate


#### get_xy_plane_atom_3d_array(one_atom)
Input: Takes in one atom feature matrix.
Takes first three features for one atom and gives xyz coordinates of the two atoms used to define the x-axis, and the
xy plane respectively


#### stack_one_atom_xyz_3D_arrays()
Iterates over all the atoms in the molecule. Every atom can be used as center for ALF. 
Stacks together 3D array for xy plane atoms, as well as 3D array that defines rest of atoms 
in xyz coordinates.

Returns all_atom_4D_array , a 4D array of shape (N_atoms, N_atoms-1, N_points, 3)
If an atom is the ALF center, all other atoms have to be represented as xyz coordinates (thus the N_atoms-1)
Every point for all other atoms has x, y, and z coordinates (so this is there the 3 comes in)

Example: Methanol has 6 atoms, when an atom is used as the ALF center only 5 atoms remain to be expressed
as xyz coordinates. These atoms have n number of points (number of timesteps in trajectory .xyz file), with
each point having x,y,z coordinates.


### ALFVisualizer.features_and_atom_names(xyz_file: str) -> (<class 'numpy.ndarray'>, <class 'list'>, <class 'list'>)
Returns features as 3D array, [atom][point][feature]
Example: 10 points water xyz file would have shape (3, 10, 3) where 3 is the number of atoms,
10 is the number of points, and 3 is the number of features


* **Parameters**

    **xyz_file** – An xyz file format file that contain trajectory information



### ALFVisualizer.resource_path(relative_path)
Get absolute path to resource, works for dev and for PyInstaller
see [https://stackoverflow.com/a/37920111](https://stackoverflow.com/a/37920111)
