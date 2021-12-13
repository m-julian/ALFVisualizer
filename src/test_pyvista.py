from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from qtpy import uic
from alfvis_core import Trajectory
from typing import List, Tuple, Dict, Union
import numpy as np
import string
import os
import sys
from pathlib import Path
from alfvis_core.constants import random_colors, default_atom_colors
from alfvis_core.start_alf_vis import _start_alf_vis_ui
from pyvista import _vtk

# Setting the Qt bindings for QtPy 5, change if using Pyside 2
# os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
os.environ["QT_API"] = "pyqt5"


# ##############################################################################
#                               XYZ FILE PARSING
# TAKES .XYZ TRAJECTORY FILE AND CALCULATES FEATURES OF EACH ATOM IN EACH POINT
###############################################################################

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

    errors_for_properties = trajectory.properties_error  # dictionary of errors for each property

    return features, atom_names, atom_name_priorities, atomic_local_frame_dict, errors_for_properties

##############################################################################################
#                     CREATING 4D ARRAY TO PLOT, SHAPE shape (N_atoms, N_atoms-1, N_points, 3)
#
##############################################################################################


class XYZArrays:

    """ Class for converting to Cartesian space.
    Creates a 3D array for each atom on which the ALF is centered (1 4D array total). Each 2D array in the 3D array
    consists of N_pointsx3 (because each point has x,y,z coords in 3D space) matrices, where each matrix contains
    the xyz coordinates of every atom that is NOT the atom on which the ALF is centered."""

    def __init__(self, all_atom_features: np.ndarray):

        self.n_atoms, self.n_points, self.n_features = all_atom_features.shape
        self.all_atom_4d_array = self.stack_one_atom_xyz_3D_arrays(all_atom_features)

    def stack_one_atom_xyz_3D_arrays(self, all_atom_features):
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

        for one_atom_features in all_atom_features:

            # xy_atom_3d_array, and polar_atoms_3d_array are both 3D atom arrays
            # (xy_atom_matrix is 2xN_pointsx3, and polar atom matrix
            # is of shape N_remaining_atomsxN_pointsx3 where N remaining atoms is N_atoms-3)

            xy_atom_3d_array = self.get_xy_plane_atom_3d_array(one_atom_features)

            # only add polar atoms if there are over 3 atoms in system
            if self.n_atoms > 3:
                polar_atoms_3d_array = self.get_polar_atom_3d_array(one_atom_features)
                # so now we stack these matrices into one 3D array that is the xyz coordinates
                # for all atoms OTHER than the atom on which the ALF is centered,
                # shape ((2+N_remaining_atoms),N_points,3)
                one_atom_total_array = np.concatenate((xy_atom_3d_array, polar_atoms_3d_array), axis=0)
                all_other_atom_3D_arrays.append(one_atom_total_array)
            else:
                one_atom_total_array = xy_atom_3d_array
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


class VisualizationWindow(QMainWindow):
    """ handles GUI and connects user commands with what to plot on pyvista plot
    see https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog videos for info on using Qt with python """

    def __init__(self, xyz_file):

        super().__init__()

        # errors_for_properties could be None in case no per-property data is read in from xyz comment line
        all_atom_features, atom_names, atoms_names_priorities, atomic_local_frame_dict, errors_for_properties = features_and_atom_names(xyz_file)

        # list of atom names
        self.atom_names: list = atom_names
        # a list of lists which gives the full ALF of every atom (not just the central, x, xy atoms)
        self.atoms_names_priorities: List[list] = atoms_names_priorities
        # dict of key: central atom, val: list of ALF indeces
        self.alf_dict: dict = atomic_local_frame_dict

        # dictionary of dictionaries, ordered as  eg. {"C1":{"O3":xyz_array, "H2":xyz_array, "H4":xyz_array ....}
        self.all_atom_dict = self.calculate_all_atom_dictionary(all_atom_features)

        self.n_timesteps = all_atom_features.shape[1]
        self._current_noncentral_data = self.all_atom_dict[atom_names[0]]

        # errors_for_properties is a list of dictionaries. Each dictionary contains property errors for each atom. Could be None if data is not in xyz.
        self.errors_for_properties = errors_for_properties
        self.cmap_properties = errors_for_properties[0].keys() if errors_for_properties is not None else None

        # current colors used for atoms
        self.current_atom_colors = dict(zip(atom_names, random_colors))
        # a saved list (in case the default colors was used and then reverted)
        self.saved_atom_colors = dict(zip(atom_names, random_colors))

        # used in initializing values for slider, atom selecter, and atom color parts, and grid
        self.checkboxes = None
        self.color_buttons_dict = None
        self.color_labels_dict = None

        ################################################
        # start user interface
        ################################################
        # load in ui
        ui_path = os.path.join(".", "test_ui.ui")
        Ui_MainWindow, _ = uic.loadUiType(ui_path)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        # initialize what is being plotted
        self._initialize_alf_vis_ui()

    def calculate_all_atom_dictionary(self, all_atom_features):
        """Returns a dictionary of dictonaries, with outer keys being the names of a central atom, eg. C1, H2, etc. The value corresponding to that outer key
        is another dictionary containing the xyz coordinates of every non-central atom that needs to be plotted. The values of the inner dictionary is a 2d-array
        of the positions of the non-central atoms for every timestep (shape n_timesteps x 3)

        :param all_atom_features: numpy array of shape n_atoms x n_timesteps x n_features
        :return: a dictonary of dictionaries containing the corresponding non-central data to plot for every central atom
        :rtype: dict
        """

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
        for center_atom, center_atom_xyzs, atom_names_prio in zip(self.atom_names, system_as_xyz.all_atom_4d_array, self.atoms_names_priorities):
            # C1  #C1 non central atom 3D array # O3, H2, H4, H5, H6 # 2, 1, 3, 4, 5
            for idx, atom_name_prio in enumerate(atom_names_prio):
                #  0 O3, 1 H2, etc. used to get out the individual 2d arrays which are ordered as:
                # x axis, xy plane, and then all atoms in polar coords
                xyz_dict[atom_name_prio] = center_atom_xyzs[idx]

            total_dict[center_atom] = xyz_dict
            xyz_dict = {}

        return total_dict

    @property
    def current_central_atom_name(self) -> str:
        """ returns the name of the current central atom. Need try except in case ui is not launched yet"""
        return self.ui.atom_names_combo.currentText()

    @property
    def current_noncentral_atom_names(self) -> list:
        """ Returns the non central atom names as a list"""
        return [name for name in self.atom_names if name != self.current_central_atom_name]

    @property
    def current_alf_str(self) -> str:
        return ', '.join(str(x) for x in self.alf_dict[self.current_central_atom_name])

    @property
    def current_selected_property(self) -> Union[str, None]:
        """ returns the name of the current selected property or `None` if no properties were read in from the xyz file."""
        if self.cmap_properties:
            return self.ui.properties_cmap_combo_box.currentText()

    @property
    def current_errors_list(self) -> list:

        if self.errors_for_properties:
            return [timestep[self.current_selected_property][self.current_central_atom_name] for timestep in self.errors_for_properties]
        return

    @property
    def center(self) -> np.ndarray:
        return np.array([0, 0, 0])

    @property
    def all_noncentral_data(self) -> dict:
        """A dictonary of non central data to plot for the current central atom."""
        return self.all_atom_dict[self.current_central_atom_name]

    @property
    def default_atom_colors(self) -> dict:
        """Returns a dictionary of default colors for the atoms in the system. (e.g. carbon is grey, nitrogen is blue)

        :rtype: dict
        """
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
    def current_central_atom_color(self) -> str:
        """ Returns the central atom color"""
        return self.current_atom_colors[self.current_central_atom_name]

    @property
    def renderer(self):
        """ Returns the renderer instance responsible for displaying things on screen"""
        return self.plotter.renderer

    @property
    def actors(self) -> dict:
        """ Returns a dictionary of actor names as keys and VTK objects as values (these are the things that are plotted)"""
        return self.renderer.actors

    def show_actor(self, actor_name):
        """ Show the actor to screen """
        actor = self.actors.get(actor_name)
        if actor:
            actor.SetVisibility(True)

    def hide_actor(self, actor_name):
        """ Hide the actor."""
        actor = self.actors.get(actor_name)
        if actor:
            actor.SetVisibility(False)

    def add_actor(self, actor, actor_name):
        """ Add an actor. This adds it to the self.actors dict"""
        self.plotter.add_mesh(actor, render_points_as_spheres=True, name=actor_name, reset_camera=False)

    def remove_actor(self, actor_name):
        """ Remove the actor. This deletes it from the self.actors dict"""
        self.plotter.remove_actor(actor_name)

    def remove_all_plotted_atoms(self):
        """ Clear the whole pyvista screen of any actors, but leaves grid and background."""
        for actor_name in list(self.actors):
            if actor_name in self.atom_names:
                self.remove_actor(actor_name=actor_name)

    def remove_grid(self):
        """ Remove the grid. """
        self.plotter.remove_bounds_axes()

    def show_grid(self):
        """ Add the grid actor and show it."""
        self.plotter.show_grid()

    def _initialize_alf_vis_ui(self):
        """ Initializes pyvista plot and user ui, with first atom ALF displayed
        Methods starting with _ only called here"""
        _start_alf_vis_ui(self)
        # called here to initialize the plot
        self.update_central_atom_and_plot()

    def update_central_atom_and_plot(self):
        """
        method runs ONLY when the central atom is changed in the GUI
        Updates central atom (always at 0,0,0 but can update color if different atom) as
        well as updates non central atom data.
        """
        self.remove_all_plotted_atoms()

        self.update_central_atom_data()
        self.update_checkboxes_box()
        self.update_colors_buttons_box()

        # self.ui.random_colors_radio.toggle()
        # self.ui.plot_all_points_radio.toggle()

        self.plot_updated_data()

    ####################################################################
    # COLORING
    ####################################################################

    def use_random_colors(self):
        """sets random colors for atoms"""
        from alfvis_core.useful_funcs import string_to_rgb

        # enable color area (if previous selected option was cmap)
        self.ui.atom_color_scroll_area.setEnabled(True)

        self.current_atom_colors = self.saved_atom_colors
        for atom in self.atom_names:
            atom_color = self.current_atom_colors[atom]
            rgb_color = string_to_rgb(atom_color)
            self.actors[atom].GetProperty().SetColor(rgb_color)
        for atom_name, button in self.color_buttons_dict.items():
            button.setStyleSheet(f"background-color : {self.current_atom_colors[atom_name]}; color: {self.current_atom_colors[atom_name]};")

        self.plot_updated_data()

    def use_default_colors(self):
        """ sets default colors for atoms"""
        from alfvis_core.useful_funcs import string_to_rgb

        self.ui.atom_color_scroll_area.setEnabled(True)
        self.current_atom_colors = self.default_atom_colors
        for atom in self.atom_names:
            atom_color = self.current_atom_colors[atom]
            rgb_color = string_to_rgb(atom_color)
            self.actors[atom].GetProperty().SetColor(rgb_color)
        for atom_name, button in self.color_buttons_dict.items():
            button.setStyleSheet(f"background-color : {self.default_atom_colors[atom_name]}; color: {self.default_atom_colors[atom_name]};")

        self.plot_updated_data()

    def use_cmap(self):
        """ Used to remove other checkboxes that cannot be used at the same time, as well as to plot the cmap."""
        # enable the plot all points radiobutton, because we are using `clicked`, this will not call `plot_all_points` automatically
        self.ui.plot_all_points_radio.setChecked(True)
        self.plot_all_points()
        self.ui.atom_color_scroll_area.setEnabled(False)
        # using cmap needs to make new actors and plot them, so it is not directly updating vtk objects.
        # possibly a way to apply a colormap directly to current vtk objects without making new ones, but that will be a lot of work and little reward

    def update_central_atom_data(self):
        """ method used to update the central ALF atom and the noncentral data associated with it, depending on selected atom in combo box"""
        self._current_noncentral_data = self.all_noncentral_data
        self.ui.atomic_local_frame.setText(self.current_alf_str)

    def plot_all_points(self):
        """ Plots data for all non-central atoms and disables the individual points widget."""

        self._current_noncentral_data = self.all_noncentral_data
        # disable individual points widget so individual points cannot be plotted.
        self.disable_plot_individual_points()
        # need to plot updated data now, so needs to be called
        self.plot_updated_data()
        # self.renderer.camera_position = "yz"
        self.renderer.camera.elevation = 0
        self.renderer.camera.azimuth = 0

    def disable_plot_individual_points(self):
        """ disable plot individual points"""
        self.ui.individual_points_widget.setEnabled(False)

    def enable_plot_individual_points(self):
        """ enable plot individual points """

        self.ui.default_atom_colors_radio.setChecked(True)
        self.use_default_colors()

        self.ui.individual_points_widget.setEnabled(True)
        self.update_individual_point_slider_status_and_box()

    def update_individual_point_slider_status_and_box(self):
        """
        Updates the status of the individual point slider depending on the checked state of the individual point checkbox.
        If that checkbox is enabled, then only 1 point will be plotted at a time. The slider can then be slider to the point of interest
        in the trajectory.
        """
        # update slider box value depending on slider position
        current_point = self.ui.individual_point_slider.value()
        self.ui.individual_point_box.setText(f"{current_point}")

        # only get one point in the trajectory corresponding to the timestep selected by the slider/box
        self._current_noncentral_data = {}
        for atom in self.all_noncentral_data.keys():
            self._current_noncentral_data[atom] = self.all_noncentral_data[atom][current_point]

        if self.current_errors_list:
            # get a list of integers for the current atom and property that are selected in combo boxes
            self.ui.property_value_for_current_point.setText(f"{self.current_errors_list[current_point]:.8f}")
        else:
            self.ui.property_value_for_current_point.setText("Not read in.")

        self.plot_updated_data()
        # resets camera to view across xz plane, see also self.renderer.reset_camera()
        # self.renderer.view_yz()
        # self.renderer.camera_position = "yz"
        # self.renderer.camera.azimuth = 45
        self.renderer.camera.elevation = 0
        self.renderer.camera.azimuth = 0
        self.renderer.Modified()

    def update_individual_point_slider_value_with_box(self):
        """ Updates the slider value based on the value of the individual point text box. The slider is updates, so the data is also automatically
        updated. No need to update data to be plotted."""

        val = "".join([s for s in self.ui.individual_point_box.text() if s.isdigit()])
        if not val:
            val = 0
        current_box_val = int(val)
        self.ui.individual_point_slider.setValue(current_box_val)  # self.update_data_and_plot() called here

    def grid_status(self):
        """ show or remove grid on pyvista plot depending on grid checkbox, updates atom data to plot"""
        if self.ui.show_grid_checkbox.isChecked():
            self.show_grid()
        elif not self.ui.show_grid_checkbox.isChecked():
            self.remove_grid()

    def update_checkboxes_box(self):
        """ Used to dynamically generate the non-central atom checkboxes.
        They can be used to plot individual noncentral atoms instead of all noncentral atoms."""

        # TODO: remove checkboxes here
        if self.checkboxes:
            for check in self.checkboxes:
                self.ui.gridLayout.removeWidget(check)
                check.deleteLater()
                check = None
            self.checkboxes = []

        self.checkboxes = []
        row = 0
        column = 0
        for atom_name in self.current_noncentral_atom_names:

            checkbox = QtWidgets.QCheckBox(f"{atom_name}")
            checkbox.setCheckState(QtCore.Qt.Checked)
            # connect each checkbox to be able to change noncentral data that is plotted
            checkbox.stateChanged.connect(self.show_or_hide_atom)
            self.checkboxes.append(checkbox)
            self.ui.gridLayout.addWidget(checkbox, row, column)
            column += 1
            # if there are 3 checkboxes on 1 row, go to the next one
            if column % 3 == 0:
                row += 1
                column = 0

    def update_colors_buttons_box(self):
        """ updates atoms that are in the color box, depending on which central atom is chosen and also which
        checkboxes are ticked in the checkbox widget."""

        # clear color buttons and labels if central atom is changed
        if self.color_buttons_dict:
            for button in self.color_buttons_dict.values():
                self.ui.gridLayout.removeWidget(button)
                button.deleteLater()
                button.deleteLater()
                button = None
        if self.color_labels_dict:
            for label in self.color_labels_dict.values():
                self.ui.gridLayout.removeWidget(label)
                label.deleteLater()
                label.deleteLater()
                label = None

        self.color_buttons_dict = {}
        self.color_labels_dict = {}

        row = 0
        column_button = 0
        column_label = 1

        for atom_name in self.atom_names:

            push_button = QtWidgets.QPushButton(atom_name)
            push_button.setStyleSheet(f"background-color : {self.current_atom_colors[atom_name]}; color: {self.current_atom_colors[atom_name]};")
            push_button.clicked.connect(self.change_atom_color_color_dialog)

            push_button_label = QtWidgets.QLabel(atom_name)

            self.color_buttons_dict[atom_name] = push_button
            self.color_labels_dict[atom_name] = push_button_label

            self.ui.gridLayout_2.addWidget(push_button, row, column_button)
            self.ui.gridLayout_2.addWidget(push_button_label, row, column_label)

            column_button += 2
            column_label += 2

            if column_button % 3 == 0:
                row += 1
                column_button = 0
                column_label = 1

    def show_or_hide_atom(self):
        """ Hide the atom from the plotter, as well as hide the atom color box that can be used to change atom color. This gets called
        by every checkbox, which has a different atom name, so then if the checkbox is checked, the data corresponding to the checkbox
        name is shown, otherwise the data is hidden."""

        sender_checkbox_text = self.sender().text()

        if self.sender().isChecked():
            # show actor
            self.show_actor(sender_checkbox_text)
            # show colorbox
            self.show_colorbutton(sender_checkbox_text)
        else:
            # hide actor
            self.hide_actor(sender_checkbox_text)
            # hide colorbox
            self.hide_colorbutton(sender_checkbox_text)

    def hide_all_atoms(self):
        """ Hides all atoms by recursively unchecking all checkboxes (when a checkbox is unchecked the `self.show_or_hide_atom` method is called.)"""
        for checkbox in self.checkboxes:
            checkbox.setCheckState(QtCore.Qt.Unchecked)

    def show_all_atoms(self):
        """ Shows all atoms by recursively checking all checkboxes (when a checkbox is checked the `self.show_or_hide_atom` method is called.)"""
        for checkbox in self.checkboxes:
            checkbox.setCheckState(QtCore.Qt.Checked)

    def show_colorbutton(self, sender_checkbox_text):
        self.color_buttons_dict[sender_checkbox_text].show()

    def hide_colorbutton(self, sender_checkbox_text):
        self.color_buttons_dict[sender_checkbox_text].hide()

    def change_atom_color_color_dialog(self):
        """ opens up color dialog and lets user select a new color for a particular atom.
        GUI is automatically frozen until user closes the color dialog box."""
        from alfvis_core.useful_funcs import hex_to_rgb

        # this opens up a color dialog box and returns a QColor Object
        color = QtWidgets.QColorDialog.getColor()
        # if the user has pressed cancel, the return color is not valid, so this does not run
        if QtGui.QColor.isValid(color):

            self.sender().setStyleSheet("")
            self.sender().setStyleSheet(f"background-color : {color.name()}; color: {color.name()};")
            self.current_atom_colors[f"{self.sender().text()}"] = color.name()
            self.saved_atom_colors[f"{self.sender().text()}"] = color.name()

            rbg_val = hex_to_rgb(color.name())
            # self.actors[f"{self.sender().text()}"].GetProperty().SetColor(rbg_val)
            # prop = _vtk.vtkProperty()
            # self.actors[f"{self.sender().text()}"].SetProperty(prop)
            # prop.SetRenderPointsAsSpheres(True)
            # prop.SetColor(rbg_val)
            self.actors[f"{self.sender().text()}"].GetProperty().SetColor(rbg_val)
            self.renderer.Modified()

    def plot_updated_data(self):
        """ plots all the data after all the checkboxes/sliders/colors etc. have been processed"""

        center = pv.PolyData(self.center)
        self.plotter.add_mesh(center, color=self.current_central_atom_color, point_size=32, render_points_as_spheres=True, name=self.current_central_atom_name)

        current_datablock = pv.MultiBlock(self._current_noncentral_data)

        if not self.ui.cmap_radio.isChecked():
            for block in current_datablock.keys():
                self.plotter.add_mesh(current_datablock[block], color=self.current_atom_colors[block], point_size=12, render_points_as_spheres=True, name=block)

        elif self.ui.cmap_radio.isChecked():
            prop = self.current_selected_property
            for block in current_datablock.keys():
                current_datablock[block][prop] = self.current_errors_list
                self.plotter.add_mesh(current_datablock[block], scalars=prop, cmap="jet", point_size=15, render_points_as_spheres=True, name=block)


def open_file_dialog(name="Select a file", default_folder=str(Path.cwd()), files_to_look_for="All Files (*);;XYZ files (*.xyz)"):
    xyz_file, _ = QtWidgets.QFileDialog.getOpenFileName(None, name, default_folder, files_to_look_for)
    return xyz_file


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    xyz_file = open_file_dialog(name="Select xyz file", default_folder="../xyzs")
    ex = VisualizationWindow(xyz_file)
    ex.show()
    app.exec_()
