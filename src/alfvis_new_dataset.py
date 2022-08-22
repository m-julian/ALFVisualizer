import pyvistaqt
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QWidget
import pyvista as pv
from qtpy import uic
from alfvis_core.trajectory import ALFVisTrajectory
from typing import List, Tuple, Dict, Union
import numpy as np
import string
from alfvis_core.constants import random_colors, default_atom_colors
from alfvis_core.start_alf_vis import _start_alf_vis_ui
from ichor.core.calculators import calculate_alf_atom_sequence, calculate_alf_features, alf_features_to_coordinates
from ichor.core.atoms import ALF
from pathlib import Path
from uuid import uuid4
import random

#########################################################################
#                       PYVISTA/ QT PLOTTING TOOL
#########################################################################


class DatasetWidget(QWidget):
    """ handles GUI and connects user commands with what to plot on pyvista plot
    see https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog videos for info on using Qt with python """

    def __init__(self, xyz_file, plotter: pyvistaqt.QtInteractor):

        super().__init__()

        # store xyz file from which features can be calculated
        self.trajectory = ALFVisTrajectory(xyz_file)
        # initialize alf
        self.initial_alf_dict = self.trajectory.alf_dict(calculate_alf_atom_sequence)
        # list of atom names
        self.atom_names: list = self.trajectory.atom_names
        # errors_for_properties is a list of dictionaries. Each dictionary contains property errors for each atom. Could be None if data is not in xyz.
        self.properties = self.trajectory.properties
        self.cmap_properties = (
            self.properties.keys()
            if self.properties is not None
            else None
        )
        # start random colors for atoms
        self.current_atom_colors = dict(
            zip(self.atom_names, random.choices(random_colors, k=self.n_atoms))
        )

        # used in initializing values for slider, atom selecter, and atom color parts, and grid
        self.checkboxes = None
        self.color_buttons_dict = None
        self.color_labels_dict = None

        # give this widget access to the plotter, so the plotter can be modified from each instance separately
        self.plotter = plotter
        # give each instance a unique id, so that each widget has its own associated data
        self.uuid = str(uuid4())

        # load in ui
        ui_path = Path(__file__).parent.absolute() / "new_dataset.ui"
        uic.loadUi(ui_path, self)
        # initialize what is being plotted
        self._initialize_alf_vis_ui()

    @property
    def atom_names_with_uuid(self):
        return [
            self.get_atom_name_with_uuid(atom_name) for atom_name in self.atom_names
        ]

    @property
    def n_atoms(self):
        return len(self.atom_names)

    @property
    def n_timesteps(self):
        return len(self.trajectory)

    @property
    def current_central_atom_name(self) -> str:
        """ returns the name of the current central atom. Need try except in case ui is not launched yet"""
        return self.atom_names_combo.currentText()

    @property
    def current_noncentral_atom_names(self) -> list:
        """ Returns the non central atom names as a list"""
        return [
            name for name in self.atom_names if name != self.current_central_atom_name
        ]

    @property
    def current_x_axis_atom_name(self) -> str:
        """ returns the name of the current central atom. Need try except in case ui is not launched yet"""
        return self.x_axis_atom_combo_box.currentText()

    @property
    def current_xy_plane_atom_name(self) -> str:
        """ returns the name of the current central atom. Need try except in case ui is not launched yet"""
        return self.xy_plane_atom_combo_box.currentText()

    @property
    def current_selected_property(self) -> Union[str, None]:
        """ returns the name of the current selected property or `None` if no properties were read in from the xyz file."""
        if self.cmap_properties:
            return self.properties_cmap_combo_box.currentText()

    @property
    def current_properties_list(self) -> list:

        if self.properties:
            return [
                timestep[self.current_selected_property][self.current_central_atom_name]
                for timestep in self.properties
            ]

    @property
    def center(self) -> np.ndarray:
        return np.array([0, 0, 0])

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

    def all_noncentral_data(self, central_atom_alf, central_atom_features) -> dict:
        """A dictionary of non central data to plot for the current central atom."""
        return self.calculate_non_central_dictionary(central_atom_alf, central_atom_features)

    def calculate_non_central_dictionary(self, central_atom_alf: ALF, central_atom_features: np.ndarray) -> dict:
        """Returns a dictionary of key: atom name, val: np array. The np array contains the xyz coordinates of every non-central atom that needs to be plotted. The values of the inner dictionary is a 2d-array
        of the positions of the non-central atoms for every timestep (shape n_timesteps x 3)

        :param central_atom_alf: The alf that is being used for the central atom
        :param central_atom_features: numpy array of shape n_timesteps x n_features
            these are the features for the current central atom
        :param central atom name: The 
        :return: a dictonary of dictionaries containing the corresponding non-central data to plot for every central atom
        :rtype: dict
        """

        # before, ordering is 0,1,2,3,4,5,...,etc.
        # after calculating the features and converting back, the order is going to be
        # central atom idx, x-axis atom index, xy-plane atom index, rest of atom indices
        # therefore, need to change the ordering of the atoms in the new xyz array to match the atom names

        # gives 3D numpy array of shape n_timesteps x n_atoms x 3
        xyz_centered_on_central_atom = alf_features_to_coordinates(central_atom_features)
        previous_atom_ordering = list(range(self.n_atoms))
        current_atom_ordering = list(central_atom_alf) + [i for i in range(self.n_atoms) if i not in central_atom_alf]
        # this will get the index that the atom was moved to after reordering.
        reverse_alf_ordering = [
            current_atom_ordering.index(num) for num in range(self.n_atoms)
        ]
        # reverse the ordering, so that the rows are the same as before
        # can now use the atom names as they were read in in initial Trajectory/PointsDirectory instance.
        xyz_centered_on_central_atom[:, previous_atom_ordering, :] = xyz_centered_on_central_atom[
            :, reverse_alf_ordering, :
        ]

        # get array of shape n_atoms x n_timesteps x 3 instead of n_timesteps x n_atoms x 3
        xyz_centered_on_central_atom = np.swapaxes(xyz_centered_on_central_atom, 0, 1)

        # use a dictionary to store atom names and their coordinates key:atom_name, val: np array of shape n_timesteps x 3
        xyz_dict = {atom_name: atom_coords for atom_name, atom_coords in zip(self.atom_names, xyz_centered_on_central_atom)}

        # delete the dictionary for the central atom as that is always 0,0,0 for every timestep
        del xyz_dict[self.current_central_atom_name]

        return xyz_dict

    def get_alf_and_features_for_central_atom(self) -> tuple:
        current_central_atom_idx = int(self.current_central_atom_name.strip(string.ascii_letters)) - 1
        currrent_x_axis_idx = int(self.current_x_axis_atom_name.strip(string.ascii_letters))  - 1

        # if molecule has 2 atoms, then there is no xy plane atom
        currrent_xy_plane_idx = None
        if self.current_xy_plane_atom_name:
            currrent_xy_plane_idx = int(self.current_xy_plane_atom_name.strip(string.ascii_letters)) - 1

        alf = ALF(current_central_atom_idx, currrent_x_axis_idx, currrent_xy_plane_idx)
        central_atom_features = self.trajectory[self.current_central_atom_name].features(calculate_alf_features, alf)

        return alf, central_atom_features

    def get_atom_name_with_uuid(self, atom_name: str):
        """ Add the uuid to the atom name"""
        return self.uuid + "-" + atom_name

    def show_actor(self, actor_name: str):
        """ Show the actor to screen """
        actor = self.actors.get(actor_name)
        if actor:
            actor.SetVisibility(True)

    def hide_actor(self, actor_name: str):
        """ Hide the actor."""
        actor = self.actors.get(actor_name)
        if actor:
            actor.SetVisibility(False)

    def add_actor(self, actor, actor_name: str):
        """ Add an actor. This adds it to the self.actors dict"""
        self.plotter.add_mesh(
            actor, render_points_as_spheres=True, name=actor_name, reset_camera=False
        )

    def remove_actor(self, actor_name: str):
        """ Remove the actor. This deletes it from the self.actors dict"""
        self.plotter.remove_actor(actor_name)

    def remove_all_plotted_atoms(self):
        """ Clear the whole pyvista screen of any actors, but leaves grid and background."""
        for actor_name in list(self.actors):
            if actor_name in self.atom_names_with_uuid:
                self.remove_actor(actor_name=actor_name)

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

        self.plot_updated_data()

    def update_property_and_plot(self):
        """ Plots updated data if property is changed."""
        self.plot_updated_data()

    def update_atom_name_combo_boxes(self):

        self.x_axis_atom_combo_box.clear()
        self.xy_plane_atom_combo_box.clear()
        self.x_axis_atom_combo_box.addItems([name for name in self.atom_names if name != self.current_central_atom_name])
        self.xy_plane_atom_combo_box.addItems([name for name in self.atom_names if name != self.current_central_atom_name and name != self.current_x_axis_atom_name])

    def update_xy_plane_atom_name_combo_boxe(self):

        self.xy_plane_atom_combo_box.clear()
        self.xy_plane_atom_combo_box.addItems([name for name in self.atom_names if name != self.current_central_atom_name and name != self.current_x_axis_atom_name])

    ####################################################################
    # COLORING
    ####################################################################

    def use_random_colors(self):
        """sets random colors for atoms and plots updated data"""

        # enable individual points if cmap was used previously
        self.plot_individual_point_radio.setEnabled(True)

        # enable color area (if previous selected option was cmap)
        self.atom_color_scroll_area.setEnabled(True)

        self.current_atom_colors = dict(
            zip(self.atom_names, random.choices(random_colors, k=self.n_atoms))
        )

        for atom_name, button in self.color_buttons_dict.items():
            button.setStyleSheet(
                f"background-color : {self.current_atom_colors[atom_name]}; color: {self.current_atom_colors[atom_name]};"
            )

        self.plot_updated_data()

    def use_default_colors(self):
        """ sets default colors for atoms and plots updated data"""

        # enable individual points if cmap was used previously
        self.plot_individual_point_radio.setEnabled(True)

        # enable color area (if previous selected option was cmap)
        self.atom_color_scroll_area.setEnabled(True)
        self.current_atom_colors = self.default_atom_colors

        for atom_name, button in self.color_buttons_dict.items():
            button.setStyleSheet(
                f"background-color : {self.default_atom_colors[atom_name]}; color: {self.default_atom_colors[atom_name]};"
            )

        self.plot_updated_data()

    def use_cmap(self):
        """ Plots all points (because we plot cmap for all points) and enables cmap with cmap bar."""
        # enable the plot all points radiobutton, because we are using `clicked`, this will not call `plot_all_points` automatically
        self.plot_all_points_radio.setChecked(True)
        # disable the individual points as they cannot be plotted when using cmap
        self.plot_individual_point_radio.setEnabled(False)
        # plot_all_points calls plot_updated_data
        self.plot_all_points()
        self.atom_color_scroll_area.setEnabled(False)
        # using cmap needs to make new actors and plot them, so it is not directly updating vtk objects.
        # possibly a way to apply a colormap directly to current vtk objects without making new ones, but that will be a lot of work and little reward

    def central_atom_color_for_all(self):

        # enable the atom color area if that was disabled previously with cmap
        self.atom_color_scroll_area.setEnabled(True)
        # enable individual points if cmap was used previously
        self.plot_individual_point_radio.setEnabled(True)

        for atom_name in list(self.current_atom_colors.keys()):
            self.current_atom_colors[atom_name] = self.current_central_atom_color

        for atom_name, button in self.color_buttons_dict.items():
            button.setStyleSheet(
                f"background-color : {self.current_atom_colors[atom_name]}; color: {self.current_atom_colors[atom_name]};"
            )

        self.plot_updated_data()

    def update_central_atom_data(self):
        """ Used to update the central ALF atom and the noncentral data associated with it, depending on selected atom in combo box"""

        alf, feats = self.get_alf_and_features_for_central_atom()
        self._current_noncentral_data = self.all_noncentral_data(alf, feats)

    def plot_all_points(self):
        """ Plots data for all non-central atoms and disables the individual points widget."""

        alf, feats = self.get_alf_and_features_for_central_atom()

        self._current_noncentral_data = self.all_noncentral_data(alf, feats)
        # disable individual points widget so individual points cannot be plotted.
        self.disable_plot_individual_points()
        # need to plot updated data now, so needs to be called
        self.plot_updated_data()
        # self.renderer.camera_position = "yz"
        self.renderer.camera.elevation = 0
        self.renderer.camera.azimuth = 0

    def disable_plot_individual_points(self):
        """ Disable plot individual points"""
        self.individual_points_widget.setEnabled(False)

    def enable_plot_individual_points(self):
        """ Enable plot individual points widget and also sets colors to default"""

        # `checked` method is used to connect to slot, but `checked` does not call the connected method automatically (since it requires user input to call it)
        # self.default_atom_colors_radio.setChecked(True)
        # self.use_default_colors()
        # enable individual point widget and plot individual point data
        self.individual_points_widget.setEnabled(True)
        self.update_individual_point_slider_status_and_box()

    def update_individual_point_slider_status_and_box(self):
        """
        Updates the status of the individual point slider depending on the checked state of the individual point checkbox.
        If that checkbox is enabled, then only 1 point will be plotted at a time. The slider can then be slider to the point of interest
        in the trajectory.
        """
        # update slider box value depending on slider position
        current_point = self.individual_point_slider.value()
        self.individual_point_box.setText(f"{current_point}")

        # only get one point in the trajectory corresponding to the timestep selected by the slider/box
        self._current_noncentral_data = {}

        alf, feats = self.get_alf_and_features_for_central_atom()
        all_noncentral_data = self.all_noncentral_data(alf, feats)

        for atom in all_noncentral_data.keys():
            self._current_noncentral_data[atom] = all_noncentral_data[atom][current_point]

        if self.current_properties_list:
            # get a list of integers for the current atom and property that are selected in combo boxes
            self.property_value_for_current_point.setText(
                f"{self.current_properties_list[current_point]:.8f}"
            )
        else:
            self.property_value_for_current_point.setText("Not read in.")

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

        val = "".join([s for s in self.individual_point_box.text() if s.isdigit()])
        if not val:
            val = 0
        current_box_val = int(val)
        self.individual_point_slider.setValue(
            current_box_val
        )  # self.update_data_and_plot() called here

    def update_checkboxes_box(self):
        """ Used to dynamically generate the non-central atom checkboxes.
        They can be used to plot individual noncentral atoms instead of all noncentral atoms."""

        # TODO: remove checkboxes here
        if self.checkboxes:
            for check in self.checkboxes:
                self.atoms_to_plot_grid.layout().removeWidget(check)
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
            self.atoms_to_plot_grid.layout().addWidget(checkbox, row, column)
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
                self.atom_color_scroll_area_grid.layout().removeWidget(button)
                button.deleteLater()
                button = None
        if self.color_labels_dict:
            for label in self.color_labels_dict.values():
                self.atom_color_scroll_area_grid.layout().removeWidget(label)
                label.deleteLater()
                label = None

        self.color_buttons_dict = {}
        self.color_labels_dict = {}

        row = 0
        column_button = 0
        column_label = 1

        for atom_name in self.atom_names:

            push_button = QtWidgets.QPushButton(atom_name)
            push_button.setStyleSheet(
                f"background-color : {self.current_atom_colors[atom_name]}; color: {self.current_atom_colors[atom_name]};"
            )
            push_button.clicked.connect(self.change_atom_color_color_dialog)

            push_button_label = QtWidgets.QLabel(atom_name)

            self.color_buttons_dict[atom_name] = push_button
            self.color_labels_dict[atom_name] = push_button_label

            self.atom_color_scroll_area_grid.layout().addWidget(
                push_button, row, column_button
            )
            self.atom_color_scroll_area_grid.layout().addWidget(
                push_button_label, row, column_label
            )

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
        text_with_uuid = self.get_atom_name_with_uuid(sender_checkbox_text)

        if self.sender().isChecked():
            # show actor
            self.show_actor(text_with_uuid)
            # show colorbox
            self.show_colorbutton(sender_checkbox_text)
        else:
            # hide actor
            self.hide_actor(text_with_uuid)
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
            self.sender().setStyleSheet(
                f"background-color : {color.name()}; color: {color.name()};"
            )
            self.current_atom_colors[f"{self.sender().text()}"] = color.name()

            # can call self.plot_updated data here instead as well
            rbg_val = hex_to_rgb(color.name())
            self.actors[
                self.get_atom_name_with_uuid(self.sender().text())
            ].GetProperty().SetColor(rbg_val)
            self.renderer.Modified()

    def plot_updated_data(self):
        """ plots all the data after all the checkboxes/sliders/colors etc. have been processed"""

        center = pv.PolyData(self.center)
        central_atom_name_with_uuid = self.get_atom_name_with_uuid(
            self.current_central_atom_name
        )
        self.plotter.add_mesh(
            center,
            color=self.current_central_atom_color,
            point_size=32,
            render_points_as_spheres=True,
            name=central_atom_name_with_uuid,
        )

        current_datablock = pv.MultiBlock(self._current_noncentral_data)

        # default point size
        point_size = 12

        if not self.cmap_radio.isChecked():

            # make points larger for individual points to see easier against other datasets
            if self.plot_individual_point_radio.isChecked():
                point_size = 20

            for block in current_datablock.keys():
                # append uuid to atom name, so that each tab has unique actor names
                block_with_uuid = self.get_atom_name_with_uuid(block)
                self.plotter.add_mesh(
                    current_datablock[block],
                    color=self.current_atom_colors[block],
                    point_size=point_size,
                    render_points_as_spheres=True,
                    name=block_with_uuid,
                )

        elif self.cmap_radio.isChecked():
            prop = self.current_selected_property
            for block in current_datablock.keys():
                block_with_uuid = self.get_atom_name_with_uuid(block)
                current_datablock[block][prop] = self.current_properties_list
                self.plotter.add_mesh(
                    current_datablock[block],
                    scalars=prop,
                    cmap="jet",
                    point_size=point_size,
                    render_points_as_spheres=True,
                    name=block_with_uuid,
                )
