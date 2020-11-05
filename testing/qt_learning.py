import sys
import random
from qtpy import QtCore, QtWidgets, QtGui

class Widget(QtWidgets.QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        hbox = QtWidgets.QHBoxLayout(self)

        splitter1 = QtWidgets.QSplitter(self)
        splitter1.setOrientation(QtCore.Qt.Horizontal)

        left = QtWidgets.QFrame(splitter1)
        left.setFrameShape(QtWidgets.QFrame.StyledPanel)

        center = QtWidgets.QFrame(splitter1)
        center.setFrameShape(QtWidgets.QFrame.StyledPanel)

        splitter2 = QtWidgets.QSplitter(splitter1)
        sizePolicy = splitter2.sizePolicy()
        sizePolicy.setHorizontalStretch(1)

        splitter2.setSizePolicy(sizePolicy)
        splitter2.setOrientation(QtCore.Qt.Horizontal)

        top_right = QtWidgets.QFrame(splitter2)
        top_right.setFrameShape(QtWidgets.QFrame.StyledPanel)
        bottom_right = QtWidgets.QFrame(splitter2)
        bottom_right.setFrameShape(QtWidgets.QFrame.StyledPanel)

        hbox.addWidget(splitter1)
        hbox.addWidget(splitter2)

        self.setGeometry(500, 500, 750, 750)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("fusion")
    w = Widget()
    w.show()
    sys.exit(app.exec_())