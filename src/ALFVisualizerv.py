from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy import uic
from trajectory import Trajectory
from typing import List, Tuple, Dict
import numpy as np
import string
import os
import sys

from copy import copy

# Setting the Qt bindings for QtPy 5, change if using Pyside 2
# os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
os.environ["QT_API"] = "pyqt5"


# ##############################################################################
#                               XYZ FILE PARSING
# TAKES .XYZ TRAJECTORY FILE AND CALCULATES FEATURES OF EACH ATOM IN EACH POINT
###############################################################################

type2mass = {'H': 1.007825, 'He': 4.002603, 'Li': 7.016005, 'Be': 9.012182, 'B': 11.009305, 'C': 12.0,
             'N': 14.003074, 'O': 15.994915, 'F': 18.998403, 'Ne': 19.99244, 'Na': 22.989769, 'Mg': 23.985042,
             'Al': 26.981539, 'Si': 27.976927, 'P': 30.973762, 'S': 31.972071, 'Cl': 34.968853, 'Ar': 39.962383,
             'K': 38.963707, 'Ca': 39.962591, 'Sc': 44.955912, 'Ti': 47.947946, 'V': 50.94396, 'Cr': 51.940508,
             'Mn': 54.938045, 'Fe': 55.9349382, 'Co': 58.933195, 'Ni': 57.935343, 'Cu': 62.929598, 'Zn': 63.929142,
             'Ga': 68.925574, 'Ge': 73.921178, 'As': 74.921597, 'Se': 79.916521, 'Br': 78.918337, 'Kr': 83.911507}

type2rad = {'H': 0.37, 'He': 0.32, 'Li': 1.34, 'Be': 0.9, 'B': 0.82, 'C': 0.77, 'N': 0.74, 'O': 0.73, 'F': 0.71,
            'Ne': 0.69, 'Na': 1.54, 'Mg': 1.3, 'Al': 1.18, 'Si': 1.11, 'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'Ar': 0.97,
            'K': 1.96, 'Ca': 1.74, 'Sc': 1.44, 'Ti': 1.36, 'V': 1.25, 'Cr': 1.27, 'Mn': 1.39, 'Fe': 1.25,
            'Co': 1.26, 'Ni': 1.21, 'Cu': 1.38, 'Zn': 1.31, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.19, 'Se': 1.16,
            'Br': 1.14, 'Kr': 1.1}

dlpoly_weights = {"H": 1.007975, "He": 4.002602, "Li": 6.9675, "Be": 9.0121831, "B": 10.8135, "C": 12.0106,
                  "N": 14.006855, "O": 15.9994, "F": 18.99840316, "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.3055,
                  "Al": 26.9815385, "Si": 28.085, "P": 30.973762, "S": 32.0675, "Cl": 35.4515, "Ar": 39.948,
                  "K": 39.0983, "Ca": 40.078, "Sc": 44.955908, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
                  "Mn": 54.938044, "Fe": 55.845, "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
                  "Ga": 69.723, "Ge": 72.63, "As": 74.921595, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
                  "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584, "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95}


def features_and_atom_names(xyz_file: str) -> Tuple[np.ndarray, List, List, Dict]:
    """ Returns features as 3D array, [atom][point][feature]
    Example: 10 points water xyz file would have shape (3, 10, 3) where 3 is the number of atoms,
    10 is the number of points, and 3 is the number of features

    :param xyz_file: An xyz file format file that contain trajectory information

    """

    trajectory = Trajectory(xyz_file)
    features = trajectory.features  # features[point][atom][feature]
    features = np.swapaxes(features, 0, 1)  # features[atom][point][feature]

    atom_names = trajectory[0].atom_names
    numbered_priorities = trajectory.priorities
    # map the numbered priorities to actual names of atoms
    atom_name_priorities = [list(map(lambda i: atom_names[i], prio)) for prio in numbered_priorities]
    for atom_name_prio in atom_name_priorities:
        del atom_name_prio[0]  # removing central atom which is not going to be present in the 3D xyz array used in plotting
        # only non-central atoms are in this array in order x_axis atom, xyz_plane atom, followed by atoms that were in spherical polar coords

    # Indeces of this ALF start from 1 (as in the actual atom names, i.e. C1, H2, etc.). It DOES NOT start at 0.
    atomic_local_frame_dict = dict(zip(atom_names, trajectory.alf_index.tolist()))

    return features, atom_names, atom_name_priorities, atomic_local_frame_dict

########################################################################
#                     CREATING 4D ARRAY TO PLOT
#
########################################################################


class XYZArrays:
    """ Class for converting to Cartesian space.
    Creates a 3D array for each atom on which the ALF is centered (1 4D array total). Each 2D array in the 3D array
    consists of N_pointsx3 (because each point has x,y,z coords in 3D space) matrices, where each matrix contains
    the xyz coordinates of every atom that is NOT the atom on which the ALF is centered."""

    def __init__(self, all_atom_features):

        self.all_atom_features = all_atom_features
        self.n_atoms, self.n_points, self.n_features = all_atom_features.shape
        self.all_atom_4d_array = self.stack_one_atom_xyz_3D_arrays()

    def stack_one_atom_xyz_3D_arrays(self):
        """Iterates over all the atoms in the molecule. Every atom can be used as center for ALF.
        Stacks together 3D array for xy plane atoms, as well as 3D array that defines rest of atoms
        in xyz coordinates.

        Returns all_atom_4D_array , a 4D array of shape (N_atoms, N_atoms-1, N_points, 3)
        If an atom is the ALF center, all other atoms have to be represented as xyz coordinates (thus the N_atoms-1)
        Every point for all other atoms has x, y, and z coordinates (so this is there the 3 comes in)

        Example: Methanol has 6 atoms, when an atom is used as the ALF center only 5 atoms remain to be expressed
        as xyz coordinates. These atoms have n number of points (number of timesteps in trajectory .xyz file), with
        each point having x,y,z coordinates."""

        all_other_atom_3D_arrays = []

        for one_atom_features in self.all_atom_features:

            # xy_atom_3d_array, and polar_atoms_3d_array are both 3D atom arrays
            # (xy_atom_matrix is 2xN_pointsx3, and polar atom matrix
            # is of shape N_remaining_atomsxN_pointsx3 where N remaining atoms is N_atoms-3)

            xy_atom_3d_array = self.get_xy_plane_atom_3d_array(one_atom_features)
            polar_atoms_3d_array = self.get_polar_atom_3d_array(one_atom_features)

            # so now we stack these matrices into one 3D array that is the xyz coordinates
            # for all atoms OTHER than the atom on which the ALF is centered,
            # shape ((2+N_remaining_atoms),N_points,3)
            one_atom_total_array = np.concatenate((xy_atom_3d_array, polar_atoms_3d_array), axis=0)
            all_other_atom_3D_arrays.append(one_atom_total_array)

        # finally we can stack these 3D arrays into one 4D array, which will contain all the info needed
        # for plotting every atom as the ALF center. This 4D array will be stored and then if can be used
        # to quickly remove/add atoms in the visualization
        all_atom_4d_array = np.stack([i for i in all_other_atom_3D_arrays])

        return all_atom_4d_array

    def get_xy_plane_atom_3d_array(self, one_atom):

        """ Input: Takes in one atom feature matrix.
        Takes first three features for one atom and gives xyz coordinates of the two atoms used to define the x-axis, and the
        xy plane respectively
        """

        # gives an n_pointx1 vector of bond 1 lengths, need to convert to matrix with xyz coordinates (n_pointsx3)
        # y and z dimension are always 0s
        tmp_bond1 = one_atom[:, [0]]
        z = np.zeros((tmp_bond1.shape[0], 2), dtype=tmp_bond1.dtype)
        bond1 = np.concatenate((tmp_bond1, z), axis=1)

        tmp_bond2 = one_atom[:, [1]]  # bond 2 length from features (n_pointsx1)
        angle12 = one_atom[:, [2]]  # angle between bond1 and bond2 (n_pontsx1)
        x_bond2 = np.multiply(np.cos(angle12), tmp_bond2)
        y_bond2 = np.multiply(np.sin(angle12), tmp_bond2)
        # z direction is always 0
        z_bond2 = np.zeros((tmp_bond2.shape[0], 1), dtype=tmp_bond2.dtype)
        bond2 = np.concatenate((x_bond2, y_bond2, z_bond2), axis=1)

        bond1_bond2_matrix = np.stack((bond1, bond2), axis=0)

        return bond1_bond2_matrix

    def get_polar_atom_3d_array(self, one_atom):

        """ Input: Takes in one atom feature matrix.
        Every three features (after the first three features that define the xy plane) have a radius, theta, and
        phi component that defines where an atom is in xyz coordinates
        This method will return a 3D array of shape n_remaining_atoms, n_points, 3, wher 3 is because x,y,z coordinate
        """

        # firist three features account for 3 atoms (central atom, atom that defines x axis, and atom that defines
        # xy plane. Therefore do not have to iterate over them)
        n_remaining_atoms = self.n_atoms - 3
        i, j, k = 3, 4, 5
        xyz_atoms_list = []

        for _ in range(n_remaining_atoms):

            r_data = one_atom[:, [i]]  # vector of r distances (n_pointsx1)
            theta_data = one_atom[:, [j]]  # vector of thetas (n_pointsx1)
            phi_data = one_atom[:, [k]]  # vector of phis (n_pointsx1)
            xx = np.multiply(np.multiply(r_data, np.sin(theta_data)), np.cos(phi_data))
            yy = np.multiply(np.multiply(r_data, np.sin(theta_data)), np.sin(phi_data))
            zz = np.multiply(r_data, np.cos(theta_data))
            xyz_atom = np.concatenate((xx, yy, zz), axis=1)  # an N_pointsx3 matrix (storing xyz info for one atom)
            # polar_atoms_xyz = np.dstack((polar_atoms_xyz, xyz_atom), axis=2)
            xyz_atoms_list.append(xyz_atom)

            i += 3
            j += 3
            k += 3

        polar_atoms_xyz_matrix = np.stack([i for i in xyz_atoms_list])

        return polar_atoms_xyz_matrix

#########################################################################
#                       PYVISTA/ QT PLOTTING TOOL
#########################################################################


class OpenXYZFile(QtWidgets.QWidget):
    """ Open File Dialog to select XYZ file"""

    def __init__(self):

        super().__init__()

        self.open_file_dialog()
        if self.xyz_file:
            self.open_pyvista_box(self.xyz_file)
        else:
            sys.exit()

    def open_file_dialog(self):
        self.xyz_file, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select an .xyz file", "../xyzs", "All Files (*);;XYZ files (*.xyz)")

    def open_pyvista_box(self, xyz_file):

        # all_atom_features are 3D array [atom][point][feature], shape is (n_atoms, n_points, n_features)
        all_atom_features, atom_names, atoms_names_priorities, atomic_local_frame_dict = features_and_atom_names(xyz_file)

        # these lines are useful if you want to print out information on maximum difference for every feature
        # all_atoms_max = np.amax(all_atom_features, axis=1) # gives n_atoms x n_features matrix of max feature values
        # all_atoms_min = np.amin(all_atom_features, axis=1) # gives n_atoms x n_features matrix of min feature values
        # all_atom_diff = all_atoms_max - all_atoms_min # differences between min and max matrices, check if phi angles are above pi for some atom
        # print(np.amax(all_atom_diff, axis=0)) # if the max is above pi for a phi angle (every 3rd feature), then it is cyclic

        system_as_xyz = XYZArrays(all_atom_features)  # gives 4D numpy array

        total_dict = {}  # dictionary of dictionaries, ordered as {"C1":{"O3":xyz_array, "H2":xyz_array, "H4":xyz_array ....},
        # "H2": {"C1":xyz_array, "O3":xyz_array.......}........,"H6":{"O3":xyz_array, "C1":xyz_array}
        # the ordering is {central_atom1: {x_axis_atom:xyz_array, xy_plane_atom:xyz_array, polar_atom:xyz_array ..... },
        #  central_atom2:{x_axis_atom:xyz_array, xy_plane_atom:xyz_array, polar_atom:xyz_array, ..... }....}
        # this was done to keep track of colors (which cannot really be done using np arrays)

        xyz_dict = {}  # this dict keeps inner dictionaries from the total_array such as {"O3":xyz_array, "H2":xyz_array, "H4":xyz_array ....}
        # it gets reset after every iteration of the loop to move onto the next atom center

        # iterate over central atoms, their respective 3D array of other atoms as xyz coords, as well as the priorities of these xyz atoms
        # (i.e which is x axis, which is xy plane etc.)
        for center_atom, center_atom_xyzs, atom_names_prio in zip(atom_names, system_as_xyz.all_atom_4d_array, atoms_names_priorities):
            # C1  #C1 non central atom 3D array # O3, H2, H4, H5, H6 # 2, 1, 3, 4, 5
            for idx, atom_name_prio in enumerate(atom_names_prio):
                #  0 O3, 1 H2, etc. used to get out the individual 2d arrays which are ordered as:
                # x axis, xy plane, and then all atoms in polar coords
                xyz_dict[atom_name_prio] = center_atom_xyzs[idx]

            total_dict[center_atom] = xyz_dict
            xyz_dict = {}

        self.main_window = VisualizationWindow(total_dict, atom_names, atomic_local_frame_dict)
        self.main_window.show()


#########################################################################
#                       PYVISTA/ QT PLOTTING TOOL
#########################################################################

# use these random colors and the atom_colors defined below if you want to see positions of individual atoms
# (which gets lost if the same colors are specified for each type of atom)
random_colors = ['#00FFFF', '#7FFFD4', '#FFE4C4', '#000000', '#DEB887',
                 '#5F9EA0', '#7FFF00', '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C', '#00FFFF',
                 '#00008B', '#008B8B', '#B8860B', '#A9A9A9', '#006400', '#A9A9A9', '#BDB76B', '#8B008B',
                 '#556B2F', '#FF8C00', '#9932CC', '#8B0000', '#E9967A', '#8FBC8F', '#483D8B', '#2F4F4F',
                 '#2F4F4F', '#00CED1', '#9400D3', '#FF1493', '#00BFFF', '#696969', '#696969', '#1E90FF', '#B22222',
                 '#FFFAF0', '#228B22', '#FF00FF', '#DCDCDC', '#F8F8FF', '#FFD700', '#DAA520', '#808080', '#008000', '#ADFF2F',
                 '#808080', '#F0FFF0', '#FF69B4', '#CD5C5C', '#4B0082', '#FFFFF0', '#F0E68C', '#E6E6FA', '#FFF0F5', '#7CFC00',
                 '#FFFACD', '#ADD8E6', '#F08080', '#E0FFFF', '#FAFAD2', '#D3D3D3', '#90EE90', '#D3D3D3', '#FFB6C1', '#FFA07A',
                 '#20B2AA', '#87CEFA', '#778899', '#778899', '#B0C4DE', '#FFFFE0', '#00FF00', '#32CD32', '#FAF0E6', '#FF00FF',
                 '#800000', '#66CDAA', '#0000CD', '#BA55D3', '#9370DB', '#3CB371', '#7B68EE', '#00FA9A', '#48D1CC', '#C71585',
                 '#191970', '#F5FFFA', '#FFE4E1', '#FFE4B5', '#FFDEAD', '#000080', '#FDF5E6', '#808000', '#6B8E23', '#FFA500',
                 '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98', '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F', '#FFC0CB',
                 '#DDA0DD', '#B0E0E6', '#800080', '#663399', '#FF0000', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#F4A460',
                 '#2E8B57', '#FFF5EE', '#A0522D', '#C0C0C0', '#87CEEB', '#6A5ACD', '#708090', '#708090', '#FFFAFA', '#00FF7F',
                 '#4682B4', '#D2B48C', '#008080', '#D8BFD8', '#FF6347', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFFFFF', '#F5F5F5',
                 '#FFFF00', '#9ACD32']

# use these default colors for each atom type if you do not care about looking at specific atoms
default_atom_colors = {"O": "red", "H": "white", "Cl": "green", "N": "blue", "C": "grey", "S": "yellow", "P": "orange"}


class VisualizationWindowDecorators:
    """ Decorators for UI """

    @staticmethod
    def clear_plot_use_grid(original_method):
        """ remove or show grid for pyvista """

        def wrapper(instance_reference):

            if instance_reference.use_grid is True:

                instance_reference.plotter.clear()
                original_method(instance_reference)
                instance_reference.plotter.show_grid()

            elif instance_reference.use_grid is False:

                instance_reference.plotter.clear()
                original_method(instance_reference)

        return wrapper


class VisualizationWindow(QMainWindow):
    """ handles GUI and connects user commands with what to plot on pyvista plot
    see https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog videos for info on using Qt with python """

    def __init__(self, all_atom_dict, atom_names, atomic_local_frame_dict):

        super().__init__()

        self.all_atom_dict = all_atom_dict
        self.atom_names = atom_names  # list of atom names
        self.alf_dict = atomic_local_frame_dict

        self.current_atom_colors = dict(zip(atom_names, random_colors))
        self.saved_atom_colors = dict(zip(atom_names, random_colors))

        # used to initialize UI to plot first central alf atom (based on index, ex. C1, O1, etc.)
        self.current_central_atom_name = atom_names[0]

        # keeps total noncentral data that can be plotted (self.all_noncentral_data)
        # self.current_noncentral_data is actually what is plotted.
        # This needs to be done to revert back to orginal whole dataset if slider is changed back to original position
        self.current_noncentral_data = self.all_noncentral_data

        # used in initializing values for slider, atom selecter, and atom color parts, and grid
        self.checkboxes = []
        self.current_checked_atoms = []
        self.color_buttons = []
        self.color_button_labels = []
        self.use_grid = True

        # Setup for ui
        ui_path = os.path.join(".", "ALFVisualizer.ui")
        Ui_MainWindow, _ = uic.loadUiType(ui_path)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self._start_alf_vis_ui()

    @property
    def center(self):
        return np.array([0, 0, 0])

    @property
    def all_noncentral_data(self):
        return self.all_atom_dict[self.current_central_atom_name]

    @property
    def default_atom_colors(self):
        atom_colors = {}
        for atom_name in self.atom_names:
            atom_element = atom_name.rstrip(string.digits)
            # if an atom's default color is defined in the default_colors dictionary
            if atom_element in default_atom_colors.keys():
                default_color = default_atom_colors[atom_element]
            # if the atom's color is not defined in the default_colors, set color to pink
            else:
                default_color = "pink"
            atom_colors[atom_name] = default_color

        return atom_colors

    @property
    def current_central_atom_color(self):
        return self.current_atom_colors[self.current_central_atom_name]

    @property
    def current_noncentral_atom_names(self):
        return [name for name in self.atom_names if name != self.current_central_atom_name]

    def _start_alf_vis_ui(self):
        """ Initializes pyvista plot and user ui, with first atom ALF displayed
        Methods starting with _ only called here"""
        self._start_combo_central_atom_names()
        self._start_individual_point_checkbox()
        self._start_individual_point_slider()
        self._start_individual_point_textbox()
        self._start_grid_checkbox()
        self._start_default_color_checkbox()
        self._start_remove_all_atoms_button()
        self._start_pyvista_plotter()
        # called here to initialize the plot
        self.update_central_atom_and_plot()

    @VisualizationWindowDecorators.clear_plot_use_grid
    def update_central_atom_and_plot(self):
        """
        method runs ONLY when the central atom is changed in the GUI
        Updates central atom (always at 0,0,0 but can update color if different atom) as
        well as updates non central atom data.
        """
        self.update_central_atom_data()
        self.update_checkboxes_widget()
        self.update_checked_atoms()
        self.update_atom_color_box_buttons()
        self.plot_updated_data()

    @VisualizationWindowDecorators.clear_plot_use_grid
    def update_noncentral_atoms_and_plot(self):
        """ updates noncentral atom data and colors to be plotted. Called by methods that do not change the current atom. The checkboxes do not need to be
        updated when data is changed but the central atom remains the same."""
        self.update_checked_atoms()
        self.update_atom_color_box_buttons()
        self.plot_updated_data()

    def _start_combo_central_atom_names(self):
        """method initializing atom names combo box from list of atom names"""
        self.ui.atom_names_combo.addItems(self.atom_names)
        self.ui.atom_names_combo.currentIndexChanged.connect(self.update_central_atom_and_plot)

    def _start_individual_point_checkbox(self):
        """method that initializes the state of the individual atom checkbox. If the checkbox is enabled, only one point is plotted at a time.
        If it is disabled, all points are plotted at the same time. This is useful to have when you want to focus on a specific point."""
        self.ui.plot_individual_point_checkbox.setCheckState(QtCore.Qt.Unchecked)
        self.ui.plot_individual_point_checkbox.stateChanged.connect(self.update_individual_point_slider_status_and_box)

    def _start_individual_point_slider(self):
        """method that initializes slider used to plot individual points instead of whole trajectory."""
        self.ui.individual_point_slider.setEnabled(False)  # default is set to false
        self.ui.individual_point_slider.setMinimum(0)  # 0 is the first timestep
        n_timesteps = self.all_noncentral_data[self.current_noncentral_atom_names[0]].shape[0]  # getting the number of timesteps from one of the atoms
        self.ui.individual_point_slider.setMaximum(n_timesteps - 1)  # need to index starting at 0, so subtract 1
        self.ui.individual_point_slider.valueChanged.connect(self.update_individual_point_slider_plotted_data)

    def _start_individual_point_textbox(self):
        """method that initializes textbox used to display the current individual point, it can also be used to go to a specific point."""
        self.ui.individual_point_box.setEnabled(False)
        self.ui.individual_point_box.setText(f"{self.ui.individual_point_slider.value()}")
        self.ui.individual_point_box.editingFinished.connect(self.update_individual_point_slider_value_with_box)

    def _start_grid_checkbox(self):
        """ Initialize checkbox that is used to show or remove grid"""
        self.ui.show_grid_checkbox.setCheckState(QtCore.Qt.Checked)
        self.ui.show_grid_checkbox.stateChanged.connect(self.grid_status)

    def _start_default_color_checkbox(self):
        """ Initialize checkbox that is used to make atoms default colors or random colors.
        Unchecked is random color and checked is default colors."""
        self.ui.default_atom_colors_checkbox.setCheckState(QtCore.Qt.Unchecked)
        self.ui.default_atom_colors_checkbox.stateChanged.connect(self.use_default_or_random_atom_colors)

    def _start_remove_all_atoms_button(self):
        """ button that unticks all noncentral atoms"""
        self.ui.remove_all_plotted_atoms.clicked.connect(self.untick_all_noncentral_atoms_checkboxes)

    def _start_pyvista_plotter(self):
        """ method to initialize pyvista plot"""
        self.plotter = QtInteractor(self.ui.pyvista_frame)
        self.plotter.set_background("royalblue", top="aliceblue")
        self.ui.horizontalLayout_3.addWidget(self.plotter.interactor)

    def update_central_atom_data(self):
        """ method used to update the central ALF atom and the noncentral data associated with it, depending on selected atom in combo box"""
        self.current_central_atom_name = self.ui.atom_names_combo.currentText()
        self.current_noncentral_data = self.all_noncentral_data
        current_alf_str = ', '.join(str(x) for x in self.alf_dict[self.current_central_atom_name])
        self.ui.atomic_local_frame.setText(current_alf_str)

    def update_individual_point_slider_status_and_box(self):
        """
        Updates the status of the individual point slider depending on the checked state of the individual point checkbox.
        If that checkbox is enabled, then only 1 point will be plotted at a time. The slider can then be slider to the point of interest
        in the trajectory.
        """
        # enable slider and point box, update noncentral data according to slider/box
        if self.ui.plot_individual_point_checkbox.isChecked() is True:
            self.ui.individual_point_slider.setEnabled(True)
            self.ui.individual_point_box.setEnabled(True)

            # update slider box value depending on slider position
            current_point = self.ui.individual_point_slider.value()
            self.ui.individual_point_box.setText(f"{current_point}")
            # only get one point in the trajectory corresponding to the timestep selected by the slider/box

            self.current_noncentral_data = {}
            for atom in self.all_noncentral_data.keys():
                self.current_noncentral_data[atom] = self.all_noncentral_data[atom][current_point]

            print("after for loop")
            print(self.all_noncentral_data)

            self.update_noncentral_atoms_and_plot()

        # disable slider and point box, plot all timesteps
        elif self.ui.plot_individual_point_checkbox.isChecked() is False:
            self.ui.individual_point_box.setEnabled(False)
            self.ui.individual_point_slider.setEnabled(False)

            self.current_noncentral_data = self.all_noncentral_data

            self.update_noncentral_atoms_and_plot()

    def update_individual_point_slider_plotted_data(self):
        """ Makes slices of the data depending on the position of the slider or the user input in the box next to the slider."""
        current_point = self.ui.individual_point_slider.value()
        self.ui.individual_point_box.setText(f"{current_point}")
        # only get one point in the trajectory corresponding to the timestep selected by the slider/box
        self.current_noncentral_data = {}
        for atom in self.all_noncentral_data.keys():
            self.current_noncentral_data[atom] = self.all_noncentral_data[atom][current_point]

        self.update_noncentral_atoms_and_plot()

    def update_individual_point_slider_value_with_box(self):
        """ Updates the slider value based on the value of the individual point text box. The slider is updates, so the data is also automatically
        updated. No need to update data to be plotted."""
        current_val = int(self.ui.individual_point_box.text())
        self.ui.individual_point_slider.setValue(current_val)

    def grid_status(self):
        """ show or remove grid on pyvista plot depending on grid checkbox, updates atom data to plot"""
        if self.ui.show_grid_checkbox.isChecked() is True:
            self.use_grid = True
        elif self.ui.show_grid_checkbox.isChecked() is False:
            self.use_grid = False

        self.update_noncentral_atoms_and_plot()

    def use_default_or_random_atom_colors(self):
        """ updates atom colors depending on default_atom_colors_checkbox. This checkbox is initialized to unchecked state, so
        random colors are used (this is initial state). If checkbox is checked, then default colors for atoms,
        i.e. oxygen:red , hydrogen:white, etc. are used"""

        if self.ui.default_atom_colors_checkbox.isChecked() is False:  # use random colors, this is initial state
            self.current_atom_colors = self.saved_atom_colors
            # self.ui.atom_color_scroll_area.setEnabled(True)

        elif self.ui.default_atom_colors_checkbox.isChecked() is True:  # use default colors
            self.current_atom_colors = self.default_atom_colors
            # self.ui.atom_color_scroll_area.setEnabled(False) # cannot update colors since using default colors

        self.update_noncentral_atoms_and_plot()

    def update_checkboxes_widget(self):
        """ Used to dynamically generate the non-central atom checkboxes.
        They can be used to plot individual noncentral atoms instead of all noncentral atoms."""

        if self.checkboxes != []:
            for check in self.checkboxes:
                self.ui.gridLayout.removeWidget(check)
                check.deleteLater()
                check = None
            self.checkboxes = []

        self.current_checked_atoms = []
        row = 0
        column = 0
        for atom in self.current_noncentral_atom_names:
            checkbox = QtWidgets.QCheckBox(f"{atom}")
            checkbox.setCheckState(QtCore.Qt.Checked)
            # connect each checkbox to be able to change noncentral data that is plotted
            checkbox.stateChanged.connect(self.update_noncentral_atoms_and_plot)
            self.checkboxes.append(checkbox)
            self.current_checked_atoms.append(checkbox.text())
            self.ui.gridLayout.addWidget(checkbox, row, column)
            column += 1
            # if there are 3 checkboxes on 1 row, go to the next one
            if column % 3 == 0:
                row += 1
                column = 0

    def untick_all_noncentral_atoms_checkboxes(self):
        """ method called after clicking the Remove All Plotted Atoms button - this unticks all noncentral atoms. Since the sates of the
        checkboxes has changed, this will cause the plot to be updated."""
        self.current_checked_atoms = []
        for checkbox in self.checkboxes:
            checkbox.setCheckState(QtCore.Qt.Unchecked)

    def update_checked_atoms(self):
        """ Method that keeps track of which checkboxes for noncentral atoms are checked. ONLY atoms that are checked are plotted"""

        for checkbox in self.checkboxes:
            if checkbox.isChecked() is False and checkbox.text() in self.current_checked_atoms:
                self.current_checked_atoms.remove(checkbox.text())
            elif checkbox.isChecked() is True and checkbox.text() not in self.current_checked_atoms:
                self.current_checked_atoms.append(checkbox.text())

    def update_atom_color_box_buttons(self):
        """ updates atoms that are in the color box, depending on which central atom is chosen and also which
        checkboxes are ticked in the checkbox widget."""

        # clear color buttons and labels if checkboxes are changed
        if self.color_buttons != []:
            for button, button_label in zip(self.color_buttons, self.color_button_labels):
                self.ui.gridLayout.removeWidget(button)
                self.ui.gridLayout.removeWidget(button_label)
                button.deleteLater()
                button_label.deleteLater()
                button = None
                button_label = None
            self.color_buttons = []
            self.color_button_labels = []

        row = 0
        column_button = 0
        column_label = 1

        # add color button for central atom as well
        push_button = QtWidgets.QPushButton(self.current_central_atom_name)
        push_button.setStyleSheet(f"background-color : {self.current_atom_colors[self.current_central_atom_name]}; color: {self.current_atom_colors[self.current_central_atom_name]};")
        push_button.clicked.connect(self.change_atom_color_color_dialog)
        push_button_label = QtWidgets.QLabel(self.current_central_atom_name)
        self.color_buttons.append(push_button)
        self.color_button_labels.append(push_button_label)
        self.ui.gridLayout_2.addWidget(push_button, row, column_button)
        self.ui.gridLayout_2.addWidget(push_button_label, row, column_label)

        column_button += 2
        column_label += 2

        for checkbox in self.checkboxes:
            if checkbox.isChecked() is True:

                push_button = QtWidgets.QPushButton(f"{checkbox.text()}")
                push_button.setStyleSheet(f"background-color : {self.current_atom_colors[checkbox.text()]}; color: {self.current_atom_colors[checkbox.text()]};")
                push_button.clicked.connect(self.change_atom_color_color_dialog)

                push_button_label = QtWidgets.QLabel(f"{checkbox.text()}")

                self.color_buttons.append(push_button)
                self.color_button_labels.append(push_button_label)

                self.ui.gridLayout_2.addWidget(push_button, row, column_button)
                self.ui.gridLayout_2.addWidget(push_button_label, row, column_label)

                column_button += 2
                column_label += 2

                if column_button % 3 == 0:
                    row += 1
                    column_button = 0
                    column_label = 1

    def change_atom_color_color_dialog(self):
        """ opens up color dialog and lets user select a new color for a particular atom.
        GUI is automatically frozen until user closes the color dialog box."""

        color = QtWidgets.QColorDialog.getColor()  # this opens up a color dialog box and returns a QColor Object
        # if the user has pressed cancel, this is set to False and does not run
        if QtGui.QColor.isValid(color) is True:

            # find which button called this method and change its color as needed
            for push_button in self.color_buttons:
                if push_button == QtCore.QObject.sender(self):

                    push_button.setStyleSheet("")
                    push_button.setStyleSheet(f"background-color : {color.name()}; color: {color.name()};")
                    self.current_atom_colors[f"{push_button.text()}"] = color.name()
                    self.saved_atom_colors[f"{push_button.text()}"] = color.name()

            self.update_noncentral_atoms_and_plot()

    def plot_updated_data(self):
        """ plots all the data after all the checkboxes/sliders/colors etc. have been processed"""

        center = pv.PolyData(self.center)
        self.plotter.add_mesh(center, color=self.current_central_atom_color, point_size=32, render_points_as_spheres=True)

        self.current_datablock = pv.MultiBlock(self.current_noncentral_data)

        for block in self.current_datablock.keys():
            if block in self.current_checked_atoms:
                self.plotter.add_mesh(self.current_datablock[block], color=self.current_atom_colors[block], point_size=15, render_points_as_spheres=True)

    # def plot_data_with_cmap(self):

    #     """ plots all the data after all the checkboxes/sliders/colors etc. have been processed"""

    #     center = pv.PolyData(self.center)
    #     self.plotter.add_mesh(center, color=self.current_central_atom_color, point_size=30, render_points_as_spheres=True)

    #     for idx, block in enumerate(self.current_datablock.keys(), start=2):
    #         if block in self.current_checked_atoms:
    #             # color = self.atom_colors.get(block)

    #             # self.current_datablock[block]["values"] = np.random.randn(20001)
    #             # add this to bottom line to add colors scalars="values", cmap="jet",
    #             self.plotter.add_mesh(self.current_datablock[block], point_size=10, render_points_as_spheres=True)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    ex = OpenXYZFile()
    app.exec_()
