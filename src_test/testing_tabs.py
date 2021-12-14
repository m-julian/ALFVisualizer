from PyQt5.QtWidgets import QWidget
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


class Button(QtWidgets.QPushButton):

    def __init__(self, name):

        super().__init__()
        self.setText(name)


class VisualizationWindow(QMainWindow):
    """ handles GUI and connects user commands with what to plot on pyvista plot
    see https://www.youtube.com/channel/UCj7i-mmOjLV17YTPIrCPkog videos for info on using Qt with python """

    def __init__(self):

        super().__init__()

        # Setup for ui
        ui_path = os.path.join(".", "testing_tabs.ui")
        uic.loadUi(ui_path, self)

        # if the + button is pressed, then launch a new tab
        self.tabWidget.setCurrentIndex(0)
        self.tabWidget.tabBarClicked.connect(self.tab_clicked)
        # set the current index to 0 (as there are two tabs to begin and the 1st tab is the + button)

        # show x that closes tabs
        self.tabWidget.setTabsClosable(True)
        # set a zero-sized widget on the + tab, so that the x is removed on it
        self.tabWidget.tabBar().setTabButton(1, QtWidgets.QTabBar.RightSide, QtWidgets.QWidget().resize(0, 0))
        self.tabWidget.tabBar().setTabButton(1, QtWidgets.QTabBar.LeftSide, QtWidgets.QWidget().resize(0, 0))
        self.tabWidget.tabCloseRequested.connect(self.close_tab)

    def tab_clicked(self, current_tab):

        if current_tab == self.sender().count()-1:
            self.tabWidget.insertTab(current_tab, QtWidgets.QWidget(), f"test{current_tab}")
            self.tabWidget.setCurrentIndex(current_tab)
            layout = QtWidgets.QVBoxLayout()
            layout.addWidget(Button(str(current_tab)))
            self.tabWidget.currentWidget().setLayout(layout)

    def close_tab(self, index):

        # TODO: remove atoms associated with tab from plotter before closing
        # if only one tab and the + tab are left, then close the main window
        if self.sender().count() == 2:
            self.close()
        # otherwise only close the tab
        else:
            self.tabWidget.removeTab(index)


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    ex = VisualizationWindow()
    ex.show()
    app.exec_()
