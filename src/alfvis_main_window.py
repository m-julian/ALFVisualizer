from qtpy import QtWidgets
import os
import sys
from pathlib import Path
from qtpy import uic
from pyvistaqt import QtInteractor
from alfvis_new_dataset import NewDatasetWidget

# Setting the Qt bindings for QtPy 5, change if using Pyside 2
# os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
os.environ["QT_API"] = "pyqt5"


def open_file_dialog(name="Select xyz file", default_folder=Path.cwd(), files_to_look_for="All Files (*);;XYZ files (*.xyz)"):

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
        self._first_dataset_tab_index = None
        self._last_dataset_tab_index = None

        # load in ui to be able to access attributes (open .ui file in QtDesigner for names)
        ui_path = Path.cwd() / "main_window.ui"
        uic.loadUi(ui_path.absolute(), self)

        # add the pyvista interactor into the pyvista part of the GUI
        self.plotter = QtInteractor(self.pyvista_part)
        self.plotter.set_background("royalblue", top="aliceblue")
        self.pyvista_part.layout().addWidget(self.plotter.interactor)

        # # initialize tab widget to 0th tab and connect tab clicked to slot
        # self.tabWidget.tabBarClicked.connect(self.tab_clicked)
        # self.tab_widget.setCurrentIndex(0)
        # # make tabs closable, but remove close button from the "+" tab by adding two 0-sized widgets on either side
        # self.tabWidget.setTabsClosable(True)
        # self.tabWidget.tabBar().setTabButton(0, QtWidgets.QTabBar.RightSide, QtWidgets.QWidget().resize(0, 0))
        # self.tabWidget.tabBar().setTabButton(0, QtWidgets.QTabBar.LeftSide, QtWidgets.QWidget().resize(0, 0))
        # # connect close signal to slot
        # self.tabWidget.tabCloseRequested.connect(self.close_tab)

        # show the main window
        self.show()

    def insert_new_datasets(self, datasets: list):
        idx = 0
        for dataset in datasets:
            self.insert_new_dataset(dataset, idx)
            idx += 1

    def insert_new_dataset(self, dataset_name: Path, current_tab: int):

        new_dataset_widget = NewDatasetWidget(dataset_name, self.plotter)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(new_dataset_widget)

        self.tab_widget.insertTab(current_tab, QtWidgets.QWidget(), f"test{current_tab}")
        self.tab_widget.setCurrentIndex(current_tab)
        self.tab_widget.currentWidget().setLayout(layout)

    # def tab_clicked(self, current_tab):

    #     total_tab_count = self.sender().count()

    #     # if "+" tab is clicked, insert a new tab (with new dataset)
    #     if current_tab == total_tab_count - 1:
    #         self.tabWidget.insertTab(current_tab, QtWidgets.QWidget(), f"test{current_tab}")
    #         self.tabWidget.setCurrentIndex(current_tab)
    #         layout = QtWidgets.QVBoxLayout()
    #         self.tabWidget.currentWidget().setLayout(layout)

    # def close_tab(self, index):

    #     # TODO: remove atoms associated with tab from plotter before closing

    #     # if only one tab and the + tab are left, then close the main window
    #     if self.sender().count() == 2:
    #         self.close()
    #     # otherwise only close the selected tab
    #     else:
    #         self.tabWidget.removeTab(index)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    xyz_file_list = open_file_dialog(name="Select xyz file", default_folder="../xyzs")
    # make instance of main window to be shown
    window = MainWindow()
    window.insert_new_datasets(xyz_file_list)
    # continue showing the app without blocking python thread
    app.exec_()
