# pyvista_alf
Using pyvista to visualize atomic ALF
# The main files are pyvista_features_visualization.py and pyvista_features_visualization.ui
The python code to calculate features and use pyvista is in the .py file and the actual ui is in the .ui file. The .ui file can be opened with Qt Designer.
The .ui file can also be converted to a python file if needed for some reason (currently it is loaded as the .ui file directly). This .ui file should ***not** be
changed manually, it is much easier to do it in Qt Designer.

# requirements.txt contains all of the packages needed to run and compile the script into binary

if you are not making an enviroment with the requirements.txt file packages, you can instead run

```
 conda install -c conda-forge pyvista 
```

after which

```
 conda install -c conda-forge pyvistaqt 
```

This installs all the needed packages (such as pyqt and Qt5).


# To compile the script with pyinstaller:

If pyinstaller is not yet installed:

```
 conda install -c conda-forge pyinstaller 
 ```

# Then to compile do:

 ```
pyi-makespec.exe --windowed ALFVIsualizer.py
 ```
 This will generate an ALFVisualizer.spec file which contains settings on how to compile the program. If you want to have the executable in one file,
 add the argument --onefile.

 In order for the program to compile, the .ui and other hidden libraries need to be added. This is done by
 adding these lines in the .spec file. (These lines are empty in the .spec file initially)
 For more info see https://pyinstaller.readthedocs.io/en/latest/spec-files.html
 and https://pyinstaller.readthedocs.io/en/stable/man/pyi-makespec.html

```
datas=[("ALFVisualizer.ui", ".")],
hiddenimports=['vtkmodules', 'vtkmodules.all', 'vtkmodules.util.numpy_support', 'vtkmodules.numpy_interface', 'vtkmodules.numpy_interface.dataset_adapter','vtkmodules.qt', 'vttmodules.util','vttmodules.vtkCommonCore','vttmodules.vtkCommonKitPython','vtkmodules.qt.QVTKRenderWindowInteractor'],
```

Finally save the changes to the .spec file and do

```
pyinstaller ALFVisualizer.spec
```

This should compile the script in the **dist** directory.

# If things are not working, some tips:

1. Set the debug flag to True in the .spec file (initially set to False)
2. Run the generated executable
3. Now more detailed errors on missing modules, etc. should be present in the errors after running the executable

