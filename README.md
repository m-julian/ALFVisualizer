# ALFVisualizer
Using pyvista to visualize atomic ALF

Select a .xyz file to visualize. Wait until a new window pops up with the visualization tool. It might take a while if the .xyz has a lot of timesteps.
Make sure to run the script when you are in the src directory because it uses relative paths for the .ui file.

Energies can also be read in from the comment line in the xyz file. If they are read in, a cmap checkbox can be ticked which displays a colormap of the enrgies.

# The main files are ALFVisualizer.py and ALFVisualizer.ui
The python code to calculate features and use pyvista is in the .py file and the actual ui is in the .ui file. The .ui file can be opened with Qt Designer.
The .ui file can also be converted to a python file if needed for some reason (currently it is loaded as the .ui file directly). This .ui file should ***not** be
changed manually, it is much easier to do it in Qt Designer.

Make sure the .xyz file you are reading in is in the form like

```
       4
 i =        0
  N        -2.2607932254        2.0945786163        0.0088077425
  H        -1.2164693277        2.0953397253        0.0040872488
  H        -2.6066449314        1.3328784610       -0.6163600241
  H        -2.6066479938        3.0138887811       -0.3460000529
       4
 i =        1
  N        -2.2632436218        2.0965535042        0.0083647927
  H        -1.2271833269        2.0900488087       -0.0078358832
  H        -2.5831176250        1.3176832550       -0.6042803346
  H        -2.5850039914        3.0066032814       -0.3402081673
       4
 i =        2
  N        -2.2655738644        2.0983784715        0.0079114939
  H        -1.2390034750        2.0849513037       -0.0198694922
  H        -2.5593871255        1.3043635538       -0.5917204616
  H        -2.5634208124        2.9987283195       -0.3349821297
```

If you want energies to be read in, add them to the comment line as floats. This enables the cmap checkbox.

```
       6
 i =        0        energy = 2.143
  C        -2.2784429253        0.3764716548        0.0613600732
  H        -1.1673363228        0.3433523315        0.0589783963
  O        -2.7413739998        1.2129700870       -0.9589752442
  H        -2.6620523911        0.7204890516        1.0463251752
  H        -2.6560170817       -0.6512416081       -0.1176831528
  H        -2.3804869890        2.1150657510       -0.7572121286
       6
 i =        1   energy = 2.675
  C        -2.2863781387        0.3834087333        0.0531928455
  H        -1.1942737791        0.3301522963        0.0224285724
  O        -2.7366005255        1.2110852658       -0.9560297206
  H        -2.6062423767        0.7046732375        1.0558981292
  H        -2.6462066446       -0.6636062224       -0.0697004711
  H        -2.3946286419        2.1003419571       -0.7288986971
       6
 i =        2   energy = 2.9035
  C        -2.2939466970        0.3894321047        0.0456932982
  H        -1.2146859241        0.3150121443       -0.0104391934
  O        -2.7315553621        1.2091993502       -0.9540227847
  H        -2.5510984567        0.6896104431        1.0693639362
  H        -2.6378315944       -0.6723766313       -0.0245673363
  H        -2.4098627840        2.0863041759       -0.7030934029
```

# Options
- Default Colors: Plot atoms with default colors, can still change individual atom colors using color box. Useful to have in cause you want to highlight a specific atom.
- Plot Cmap: If energies are read in from the .xyz file, this makes a cmap of all the plotted points (so can be used to plot energies/multipoles or model predictions).
- Plot Individual Point: Only 1 points is displayed at a time. Use the slider or write in box to go to specific point. Starts from 0. This cannot be used with cmap.

# Installation
if you are not making an enviroment with the requirements.txt file packages, you can instead run
```
 conda install -c conda-forge pyvista 
```

after which

```
 conda install -c conda-forge pyvistaqt 
```

This installs all the needed packages (such as pyqt and Qt5).