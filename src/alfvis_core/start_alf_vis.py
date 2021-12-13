#####################################################
# Initial Setup and Connect QObjects signals to slots
#####################################################
from pyvistaqt import QtInteractor


def _start_pyvista_plotter(self):
    """ method to initialize pyvista plot"""
    self.plotter = QtInteractor(self.ui.pyvista_frame)
    self.plotter.set_background("royalblue", top="aliceblue")
    self.ui.horizontalLayout_2.addWidget(self.plotter.interactor)


def _start_global_settings(self):
    """ Initialize checkbox that is used to show or remove grid, as well as initialize the grid (which can then be controlled by the grid checkbox).
    Any global settings which must be set for all datasets are going to be initialized here."""
    self.ui.show_grid_checkbox.stateChanged.connect(self.grid_status)
    self.grid_status()  # initialize grid here


def _start_atom_center_and_property_settings(self):
    """ Sets up the atom center and property combo boxes in the Atom Center and Property Settings GroupBox. These define what data is going to be
    available for plotting, so they must be initialized before everything else relating to plotting."""
    # initializing atom names combo box from list of atom names
    self.ui.atom_names_combo.addItems(self.atom_names)
    self.ui.atom_names_combo.currentIndexChanged.connect(self.update_central_atom_and_plot)
    # initializing the property names for which cmap can be plotted
    if self.cmap_properties:
        self.ui.properties_cmap_combo_box.addItems(self.cmap_properties)
        self.ui.properties_cmap_combo_box.currentIndexChanged.connect(self.update_property_and_plot)
    else:
        self.ui.properties_cmap_combo_box.setEnabled(False)
        self.ui.properties_cmap_combo_box.setToolTip("Energies were not read in.")


def _start_coloring_settings(self):
    """ Connect the color radio buttons to their corresponding methods. Also, initialize the colors to random

    .. note::
        `clicked` is used instead of `toggled`. `clicked` only calls the connected slot when there is user input in the gui (so if the button's check state
        is changed in code, the connected method will not run automatically). However, `clicked` is used because `toggled` causes the connected method
        to be ran even if the button is unchecked (since these are mutually exclusive buttons)
    """

    # use clicked instead of toggled because we only want to do this when the button is clicked
    self.ui.random_colors_radio.clicked.connect(self.use_random_colors)
    self.ui.default_atom_colors_radio.clicked.connect(self.use_default_colors)
    self.ui.cmap_radio.clicked.connect(self.use_cmap)


def _start_points_settings(self):
    """method that initializes the state of the individual atom checkbox. If the checkbox is enabled, only one point is plotted at a time.
    If it is disabled, all points are plotted at the same time. This is useful to have when you want to focus on a specific point.

    .. note::
        `clicked` is used instead of `toggled`. `clicked` only calls the connected slot when there is user input in the gui (so if the button's check state
        is changed in code, the connected method will not run automatically). However, `clicked` is used because `toggled` causes the connected method
        to be ran even if the button is unchecked (since these are mutually exclusive buttons)
    """

    # plot all data that has been loaded in for the current central atom
    self._current_noncentral_data = self.all_noncentral_data

    # use clicked instead of toggled because we only want to do this when the radio button is clicked
    # otherwise, because these two radio buttons are mutually exclusive, using toggle will also call the method connected to the other button (that is being deselected)
    self.ui.plot_all_points_radio.clicked.connect(self.plot_all_points)

    # set up slider values by which to index the dataset
    self.ui.individual_point_slider.setMinimum(0)  # 0 is the first timestep
    self.ui.individual_point_slider.setMaximum(self.n_timesteps - 1)  # need to index starting at 0, so subtract 1
    # set up the box that shows what index we are currently on
    self.ui.individual_point_box.setText(f"{self.ui.individual_point_slider.value()}")
    # on editing the value of this box, we can update the slider (which then updates the individual point that is plotted)
    self.ui.plot_individual_point_radio.clicked.connect(self.enable_plot_individual_points)
    self.ui.individual_point_box.editingFinished.connect(self.update_individual_point_slider_value_with_box)
    self.ui.individual_point_slider.valueChanged.connect(self.update_individual_point_slider_status_and_box)

    # by default the individual point widget is disabled (as the default radio button is the Plot All Points instead of Individual Point)
    self.ui.individual_points_widget.setEnabled(False)

    # if property data has not been read in
    if not self.cmap_properties:
        self.ui.property_value_for_current_point.setText("Not read in.")
    # if data has been read in but individual energies are still not being plotted
    else:
        self.ui.property_value_for_current_point.setText("Not available.")
        # self.ui.error_for_current_point.setText(str(self.current_errors_list[self.ui.individual_point_slider.value()]))


def _start_hide_all_atoms_button(self):
    """ button that unticks all noncentral atoms"""
    self.ui.show_all_plotted_atoms_button.clicked.connect(self.show_all_atoms)
    self.ui.hide_all_plotted_atoms_button.clicked.connect(self.hide_all_atoms)


def _start_alf_vis_ui(self):
    _start_pyvista_plotter(self)
    _start_global_settings(self)
    _start_atom_center_and_property_settings(self)
    _start_coloring_settings(self)
    _start_points_settings(self)
    _start_hide_all_atoms_button(self)
