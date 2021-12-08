import pyvista as pv
from pyvistaqt import BackgroundPlotter
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import sys

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    plotter = BackgroundPlotter(app=app)
    _ = plotter.add_mesh(pv.Sphere())
    app.exec_()
