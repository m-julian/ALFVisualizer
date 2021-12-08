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
        ui_path = os.path.join(".", "ALFVisualizer.ui")
        Ui_MainWindow, _ = uic.loadUiType(ui_path)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.plotter = QtInteractor(self.ui.pyvista_frame)
        self.plotter.set_background("royalblue", top="aliceblue")
        print(dir(self.plotter))
        self.ui.horizontalLayout_3.addWidget(self.plotter.interactor)

        circle = pv.Circle()
        self.plotter.add_mesh(circle, render_points_as_spheres=True)

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    ex = VisualizationWindow()
    ex.show()
    app.exec_()
