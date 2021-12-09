from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy import uic
import os
import sys
from pathlib import Path

# Setting the Qt bindings for QtPy 5, change if using Pyside 2
# os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
os.environ["QT_API"] = "pyqt5"


#########################################################################
#                       PYVISTA/ QT PLOTTING TOOL
#########################################################################

# use these random colors and the atom_colors defined below if you want to see positions of individual atoms
# (which gets lost if the same colors are specified for each type of atom)

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

    def __init__(self):

        super().__init__()

        # Setup for ui
        ui_path = os.path.join(".", "test.ui")
        Ui_MainWindow, _ = uic.loadUiType(ui_path)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.plotter = QtInteractor(self.ui.pyvista_frame, line_smoothing=True, point_smoothing=True)
        self.plotter.set_background("royalblue", top="aliceblue")
        print(dir(self.plotter))
        self.ui.horizontalLayout_3.addWidget(self.plotter)

        self.plotter.add_mesh(pv.Sphere(), render_points_as_spheres=True, name="sphere")
        self.plotter.add_mesh(pv.Pyramid(), render_points_as_spheres=True, name="pyramid")

        self.ui.remove_all.clicked.connect(self.remove_all_actors)
        self.ui.remove_circle.clicked.connect(self.remove_actor)
        self.ui.add_circle.clicked.connect(self.add_actor)
        self.ui.remove_grid.clicked.connect(self.remove_grid)
        self.ui.show_grid.clicked.connect(self.show_grid)
        self.ui.visibility_on_button.clicked.connect(self.visibility_on)
        self.ui.visibility_off_button.clicked.connect(self.visibility_off)

    @property
    def renderer(self):
        return self.plotter.renderer

    @property
    def actors(self):
        return self.renderer.actors

    def remove_all_actors(self):
        for actor_name in list(self.actors):
            self.remove_actor(actor_name=actor_name)

    def remove_actor(self, signal=False, actor_name="sphere"):
        self.plotter.remove_actor(actor_name)
        print("removed actor", actor_name)

    def add_actor(self):
        actor = pv.Sphere()
        actor_name = "sphere"
        self.plotter.add_mesh(actor, render_points_as_spheres=True, name=actor_name, reset_camera=False)
        print("added actor", actor, actor_name)

    def remove_grid(self):
        self.plotter.remove_bounds_axes()

    def show_grid(self):
        self.plotter.show_grid()

    def visibility_off(self):
        actor_name = "sphere"
        actor = self.actors.get(actor_name)
        if actor is not None:
            actor.SetVisibility(False)

    def visibility_on(self):
        actor_name = "sphere"
        actor = self.actors.get(actor_name)
        if actor is not None:
            actor.SetVisibility(True)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    ex = VisualizationWindow()
    ex.show()
    app.exec_()
