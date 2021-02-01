# pyvista_alf
Using pyvista to visualize atomic ALF

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


# To compile the script with pyinstaller into one file:

If pyinstaller is not yet installed:

```
 conda install -c conda-forge pyinstaller 
 ```

 Then to compile do:

 ```
 pyinstaller --onefile pyvista_features_visualization.py
 ```