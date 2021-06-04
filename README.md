# pyvista_alf
Using pyvista to visualize atomic ALF

Select a .xyz file to visualize. Wait until a new window pops up with the visualization tool. It might take a while if the .xyz has a lot of timesteps.
Make sure to run the script when you are in the src directory because it uses relative paths for the .ui file.

# The main files are ALFVisualizer.py and ALFVisualizer.ui
The python code to calculate features and use pyvista is in the .py file and the actual ui is in the .ui file. The .ui file can be opened with Qt Designer.
The .ui file can also be converted to a python file if needed for some reason (currently it is loaded as the .ui file directly). This .ui file should ***not** be
changed manually, it is much easier to do it in Qt Designer.

Make sure the .xyz file you are reading in is in the form like

```
       4
 i =        0, time =        0.000, E =       -11.6673224916
  N        -2.2607932254        2.0945786163        0.0088077425
  H        -1.2164693277        2.0953397253        0.0040872488
  H        -2.6066449314        1.3328784610       -0.6163600241
  H        -2.6066479938        3.0138887811       -0.3460000529
       4
 i =        1, time =        0.500, E =       -11.6677567782
  N        -2.2632436218        2.0965535042        0.0083647927
  H        -1.2271833269        2.0900488087       -0.0078358832
  H        -2.5831176250        1.3176832550       -0.6042803346
  H        -2.5850039914        3.0066032814       -0.3402081673
       4
 i =        2, time =        1.000, E =       -11.6671406813
  N        -2.2655738644        2.0983784715        0.0079114939
  H        -1.2390034750        2.0849513037       -0.0198694922
  H        -2.5593871255        1.3043635538       -0.5917204616
  H        -2.5634208124        2.9987283195       -0.3349821297
```


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
pyi-makespec --windowed ALFVIsualizer.py
 ```
 This will generate an ALFVisualizer.spec file which contains settings on how to compile the program. If you want to have the executable in one file,
 add the argument `--onefile `. The one file executable is around 350 mb and a bit slower to run.

 In order for the program to compile, the .ui and other hidden libraries need to be added. This is done by
 adding these lines in the .spec file. (These lines are empty in the .spec file initially)
 For more info see https://pyinstaller.readthedocs.io/en/latest/spec-files.html , https://pyinstaller.readthedocs.io/en/stable/man/pyi-makespec.html
 and https://github.com/pyvista/pyvista-support/issues/167

```
datas=[("ALFVisualizer.ui", ".")],
hiddenimports=['vtkmodules', 'vtkmodules.all', 'vtkmodules.util.numpy_support', 'vtkmodules.numpy_interface', 'vtkmodules.numpy_interface.dataset_adapter','vtkmodules.qt', 'vttmodules.util','vttmodules.vtkCommonCore','vttmodules.vtkCommonKitPython','vtkmodules.qt.QVTKRenderWindowInteractor'],
```

Finally save the changes to the .spec file and do

```
pyinstaller ALFVisualizer.spec
```

This should compile the script in the **dist** directory. If the `--onefile` option is specified, it will just show up as an executable. If it was not specified,
a folder will show up with the executable inside it. Make sure to move the whole folder if it was compiled like that.

# If things are not working, some tips:

1. Set the debug flag to True in the .spec file (initially set to False)
2. Run the generated executable
3. Now more detailed errors on missing modules, etc. should be present in the errors after running the executable

