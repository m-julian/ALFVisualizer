# Setting the Qt bindings for QtPy 5, change if using Pyside 2
#os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
import os
import sys
os.environ["QT_API"] = "pyqt5"

from qtpy import QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from pyvistaqt import QtInteractor

import numpy as np

class MainWindow(QMainWindow):

    def __init__(self, data, parent=None, show=True):

        QtWidgets.QMainWindow.__init__(self, parent)

        self.data = data

        # create the frame
        self.frame = QtWidgets.QFrame()
        vlayout = QtWidgets.QVBoxLayout()

        # add the pyvista interactor object
        self.plotter = QtInteractor(self.frame)
        vlayout.addWidget(self.plotter.interactor)

        self.frame.setLayout(vlayout)
        self.setCentralWidget(self.frame)

        # simple menu to demo functions
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        exitButton = QtWidgets.QAction('Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)

        # allow adding a sphere
        meshMenu = mainMenu.addMenu('Mesh')
        self.add_points_action = QtWidgets.QAction('Add Sphere', self)
        self.add_ponts_action2 = QtWidgets.QAction('Remove Sphere', self)
        self.add_points_action.triggered.connect(self.add_data)
        self.add_ponts_action2.triggered.connect(self.remove_data)
        meshMenu.addAction(self.add_points_action)
        meshMenu.addAction(self.add_ponts_action2)

        if show:
            self.show()

    def add_data(self):
        """ add a sphere to the pyqt frame """

        self.plotter.add_mesh(self.data, show_edges=True)
        self.plotter.reset_camera()

    def remove_data(self):
        """ add a sphere to the pyqt frame """

        self.plotter.clear()


points = np.array(([1,1,1],[2,2,2],[3,2,1]))
points2 = np.array(([7,1,1],[6,2,5],[5,2,1]))
data_cloud1 = pv.PolyData(points)
data_cloud2 = pv.PolyData(points2)

if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow(data_cloud1)
    sys.exit(app.exec_())

