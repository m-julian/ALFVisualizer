
import sys

# Setting the Qt bindings for QtPy
import os
os.environ["QT_API"] = "pyqt5"

from qtpy import QtWidgets
from qtpy.QtWidgets import QMainWindow

import numpy as np

import pyvista as pv
from pyvistaqt import QtInteractor

class MainWindow(QMainWindow):

    def __init__(self, parent=None, show=True):

        QtWidgets.QMainWindow.__init__(self, parent)

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
        self.add_sphere_action = QtWidgets.QAction('Add Sphere', self)
        self.add_sphere_action2 = QtWidgets.QAction('Remove Sphere', self)
        self.add_sphere_action.triggered.connect(self.add_sphere)
        self.add_sphere_action2.triggered.connect(self.remove_sphere)
        meshMenu.addAction(self.add_sphere_action)
        meshMenu.addAction(self.add_sphere_action2)

        if show:
            self.show()

    def add_sphere(self):
        """ add a sphere to the pyqt frame """
        sphere = pv.Sphere()
        self.plotter.add_mesh(sphere, show_edges=True)
        self.plotter.reset_camera()

    def remove_sphere(self):
        """ add a sphere to the pyqt frame """

        self.plotter.clear()


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())

