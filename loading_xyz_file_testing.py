import os
import sys
os.environ["QT_API"] = "pyqt5"
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy import uic
# import matplotlib

# app = QtWidgets.QApplication([])

# w = QtWidgets.QFileDialog()
# w.show()

# app.exec()


class App(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        self.openFileNameDialog()

    def openFileNameDialog(self):
        options = QtWidgets.QFileDialog.Options()
        options = QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Select XYZ file to visualize", "","All Files (*);;XYZ Files (*.xyz)",)  #options=options
        # if fileName:
        #     print(fileName)


app = QtWidgets.QApplication(sys.argv)
ex = App()
# ex.show()
app.exec_()

