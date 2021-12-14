from PyQt5.QtWidgets import QWidget
from qtpy import QtWidgets
import os
import sys
from pathlib import Path
from qtpy import uic
from pyvistaqt import QtInteractor
from alfvis_new_dataset import DatasetWidget
from typing import List

# Setting the Qt bindings for QtPy 5, change if using Pyside 2
# os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
os.environ["QT_API"] = "pyqt5"


def open_file_dialog(name="Select xyz file", default_folder=str(Path.cwd()), files_to_look_for="All Files (*);;XYZ files (*.xyz)"):

    # can select multiple names
    file_list, _ = QtWidgets.QFileDialog.getOpenFileNames(None, name, default_folder, files_to_look_for)

    # if there were files selected and list if not empty
    if file_list:
        xyz_file_list = [Path(f) for f in file_list if Path(f).suffix == ".xyz"]
        return xyz_file_list

    return


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()

        self._datasets = {}

        # load in ui to be able to access attributes (open .ui file in QtDesigner for names)
        ui_path = Path.cwd() / "main_window.ui"
        uic.loadUi(ui_path.absolute(), self)

        # add the pyvista interactor into the pyvista part of the GUI
        self.plotter = QtInteractor(self.pyvista_part)
        self.plotter.set_background("royalblue", top="aliceblue")
        self.pyvista_part.layout().addWidget(self.plotter.interactor)

        # initialize tab widget to 0th tab and connect tab clicked to slot
        self.tab_widget.tabBarClicked.connect(self.new_tab_clicked)
        # tabs are closable (set in .ui file), but remove close button from the "+" tab by adding two 0-sized widgets on either side
        self.tab_widget.tabBar().setTabButton(0, QtWidgets.QTabBar.RightSide, QtWidgets.QWidget().resize(0, 0))
        self.tab_widget.tabBar().setTabButton(0, QtWidgets.QTabBar.LeftSide, QtWidgets.QWidget().resize(0, 0))
        # connect close signal to slot
        self.tab_widget.tabCloseRequested.connect(self.close_tab)

        self.show_grid_checkbox.connect()


    @property
    def n_tabs(self):
        """ Returns the number of opened tabs"""
        return self.tab_widget.count()

    @property
    def last_tab_index(self):
        """ Gives the index of the last widget (which is always the `add new tab` tab). Tab indexing starts at 0."""
        return self.n_tabs - 1

    def insert_new_datasets(self, datasets: List[Path]):
        """ Inserts multiple datasets into individual new tabs. 

        :param datasets: A list of Path objects pointing to .xyz files
        """
        for dataset in datasets:
            self.insert_new_dataset(dataset, self.last_tab_index)

    def insert_new_dataset(self, dataset_path: Path, index: int):
        """ Inserts a dataset into a new tab.

        :param dataset_path: A Path to a dataset (.xyz file)
        :param index: The tab index where to insert the new dataset
        """

        # make an instance of the new dataset widget and store it, so it can be accessed by main window
        new_dataset_widget = DatasetWidget(dataset_path, self.plotter)
        # self._datasets[new_dataset_widget.uuid] = new_dataset_widget

        self.tab_widget.insertTab(index, QtWidgets.QWidget(), dataset_path.stem)
        self.tab_widget.setCurrentIndex(index)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(new_dataset_widget)
        self.tab_widget.currentWidget().setLayout(layout)

    def new_tab_clicked(self, current_tab):

        # if "+" tab is clicked, insert a new tab (with new dataset)
        if current_tab == self.last_tab_index:

            xyz_file_list = open_file_dialog()
            if xyz_file_list:
                self.insert_new_datasets(xyz_file_list)

    def close_tab(self, index: int):
        """ Closes a given tab, or if the final tab is being closed, close the application.

        :param index: The index of the tab to be closed.
        """

        # if only one tab and the + tab are left, then close the main window
        if self.sender().count() == 2:
            self.close()
        # otherwise only close the selected tab and remove the data associated with the tab
        else:
            # access the NewDatasetWidget instance for this particular tab
            widget = self.tab_widget.widget(index).findChild(DatasetWidget)
            widget.remove_all_plotted_atoms()
            self.tab_widget.removeTab(index)

    def remove_grid(self):
        """ Remove the grid. """
        self.plotter.remove_bounds_axes()

    def show_grid(self):
        """ Add the grid actor and show it."""
        self.plotter.show_grid()

    def grid_status(self):
        # TODO: move to global settings for main window
        """ show or remove grid on pyvista plot depending on grid checkbox, updates atom data to plot"""
        if self.show_grid_checkbox.isChecked():
            self.show_grid()
        elif not self.show_grid_checkbox.isChecked():
            self.remove_grid()


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    xyz_file_list = open_file_dialog()
    # make instance of main window to be shown
    main_window = MainWindow()
    if xyz_file_list:
        main_window.insert_new_datasets(xyz_file_list)
    # continue showing the app without blocking python thread
    main_window.show()
    app.exec_()
