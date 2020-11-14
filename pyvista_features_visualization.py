import re
import numpy as np
import itertools as it
import math
from functools import lru_cache
import string
import pandas as pd
from glob import glob
from itertools import cycle
from copy import copy
# Setting the Qt bindings for QtPy 5, change if using Pyside 2
#os.environ["QT_API"] = "pyside2"
# http://qtdocs.pyvista.org/usage.html
import os
import sys
os.environ["QT_API"] = "pyqt5"
from qtpy import QtCore, QtGui, QtWidgets
from qtpy.QtWidgets import QMainWindow
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy import uic
# import matplotlib

###############################################################################
#                                XYZ FILE PARSING
# TAKES .XYZ TRAJECTORY FILE AND CALCULATES FEATURES OF EACH ATOM IN EACH POINT 
###############################################################################

type2mass = {'H': 1.007825, 'He': 4.002603, 'Li': 7.016005, 'Be': 9.012182, 'B': 11.009305, 'C': 12.0,
             'N': 14.003074, 'O': 15.994915, 'F': 18.998403, 'Ne': 19.99244, 'Na': 22.989769, 'Mg': 23.985042,
             'Al': 26.981539, 'Si': 27.976927, 'P': 30.973762, 'S': 31.972071, 'Cl': 34.968853, 'Ar': 39.962383,
             'K': 38.963707, 'Ca': 39.962591, 'Sc': 44.955912, 'Ti': 47.947946, 'V': 50.94396, 'Cr': 51.940508,
             'Mn': 54.938045, 'Fe': 55.9349382, 'Co': 58.933195, 'Ni': 57.935343, 'Cu': 62.929598, 'Zn': 63.929142,
             'Ga': 68.925574, 'Ge': 73.921178, 'As': 74.921597, 'Se': 79.916521, 'Br': 78.918337, 'Kr': 83.911507}

type2rad = {'H': 0.37, 'He': 0.32, 'Li': 1.34, 'Be': 0.9, 'B': 0.82, 'C': 0.77, 'N': 0.74, 'O': 0.73, 'F': 0.71,
            'Ne': 0.69, 'Na': 1.54, 'Mg': 1.3, 'Al': 1.18, 'Si': 1.11, 'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'Ar': 0.97,
            'K': 1.96, 'Ca': 1.74, 'Sc': 1.44, 'Ti': 1.36, 'V': 1.25, 'Cr': 1.27, 'Mn': 1.39, 'Fe': 1.25,
            'Co': 1.26, 'Ni': 1.21, 'Cu': 1.38, 'Zn': 1.31, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.19, 'Se': 1.16,
            'Br': 1.14, 'Kr': 1.1}

dlpoly_weights =  {"H": 1.007975, "He": 4.002602, "Li": 6.9675, "Be": 9.0121831, "B": 10.8135, "C": 12.0106,
                   "N": 14.006855, "O": 15.9994, "F": 18.99840316, "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.3055,
                   "Al": 26.9815385, "Si": 28.085, "P": 30.973762, "S": 32.0675, "Cl": 35.4515, "Ar": 39.948,
                   "K": 39.0983, "Ca": 40.078, "Sc": 44.955908, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
                   "Mn": 54.938044, "Fe": 55.845, "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
                   "Ga": 69.723, "Ge": 72.63, "As": 74.921595, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
                   "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584, "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95}


class Trajectory:
    """ Reads trajectory file"""

    def __init__(self, fname, read=False):
        self.fname = fname
        self._trajectory = []

        if self.fname and read:
            self.read()
    
    def read(self):
        """ read in trajectory file (typically .xyz)"""
        with open(self.fname, "r") as f:
            atoms = Atoms()
            for line in f: # reads through lines in .xyz file
                if not line.strip():
                    continue # ignore blank lines
                elif re.match(r"^\s+\d+$", line): # finds number of atoms in xyz file
                    natoms = int(line)
                    while len(atoms) < natoms:
                        # add all atoms in a timestep to an instance of Atoms
                        line = next(f)
                        if re.match(r"\s*\w+(\s+[+-]?\d+.\d+([Ee]?[+-]?\d+)?){3}", line): # if re finds atom type and x, y,z coordinates it will add to list
                            atoms.add(line) # see class Atoms add method, line is a coordinate line
                    self.add(atoms)
                    atoms = Atoms()
    
    def add(self, atoms):
        """ each item in self._trajectory is of Atoms class,
        printing self._trajectory prints out based on Atoms __str__ method"""
        self._trajectory.append(atoms)

    @property
    def features(self):
        return np.array([point.features for point in self])

    @property
    def atom_names(self):
        return [atom.name for atom in self[0]]

    @property
    def nfeatures(self):
        return self[0].nfeatures

    @property
    def priorities(self):
        return(Atoms.PRIORITIES)

    def __len__(self):
        return len(self._trajectory)
    
    def __getitem__(self, i):
        return self._trajectory[i]


class Atoms:
    """ Deals with 1 trajectory timestep / all atoms in a timestep"""
    ALF = []
    PRIORITIES = [] # simmilar to ALF variable but will also store info about atoms that are NOT used as x
    # and xy plane. This will then be used to create a dictionary containing atoms and what needs to be plotted in
    # the GUI

    def __init__(self, atoms=None):
        Atom.counter = it.count(1)
        self._atoms = []

        if not atoms is None:
            self.add(atoms)

    def add(self, atom):
        """ adds an instance of Atom to _atoms"""
        if isinstance(atom, str):
            self._atoms.append(Atom(atom))
        elif isinstance(atom, Atom):
            self._atoms.append(atom)
        elif isinstance(atom, (list, Atoms)):
            for a in atom:
                self.add(a)

    @property
    def priority(self):
        return sum(self.masses)

    @property
    def max_priority(self):
        prev_priorities = []
        while True:
            priorities = [atom.priority for atom in self]
            if priorities.count(max(priorities)) == 1 or prev_priorities == priorities:
                break
            else:
                prev_priorities = priorities
        for atom in self:
            atom.reset_level()
        return self[priorities.index(max(priorities))]

    @property
    def masses(self):
        return [atom.mass for atom in self]

    @property
    def atoms(self):
        return [atom.atom_num for atom in self]

    @property
    def empty(self):
        return len(self) == 0

    def connect(self, iatom, jatom):
        iatom.set_bond(jatom)
        jatom.set_bond(iatom)

    @property
    @lru_cache()
    def connectivity(self):
        connectivity = np.zeros((len(self), len(self)))
        for i, iatom in enumerate(self):
            for j, jatom in enumerate(self):
                if not iatom == jatom:
                    max_dist = 1.2 * (iatom.radius + jatom.radius)

                    if iatom.dist(jatom) < max_dist:
                        connectivity[i][j] = 1
                        self.connect(iatom, jatom)
        return connectivity

    def to_angstroms(self):
        for atom in self:
            atom.to_angstroms()

    def to_bohr(self):
        for atom in self:
            atom.to_bohr()
    
    @property
    def alf(self):
        return [[iatom.num for iatom in atom.alf] for atom in self]

    def calculate_alf(self):
        self.connectivity
        for iatom in self:
            for alf_axis in range(2):
                queue = iatom.bonds - iatom.alf
                if queue.empty:
                    for atom in iatom.bonds:
                        queue.add(atom.bonds)
                    queue -= iatom.alf
                iatom.add_alf_atom(queue.max_priority)
        Atoms.ALF = self.alf
        Atoms.PRIORITIES = self.alf
        for alf_atom in Atoms.PRIORITIES:
            for n_atom in range(1,len(self)+1):
                if n_atom in alf_atom:
                    continue
                else:
                    alf_atom.append(n_atom)

        self.set_alf()

    def set_alf(self):
        for atom, atom_alf in zip(self, Atoms.ALF):
            atom.x_axis = self[atom_alf[1]-1]
            atom.xy_plane = self[atom_alf[2]-1]

    def calculate_features(self):
        """ calculates features for atoms in one timestep"""
        if not Atoms.ALF:
            self.calculate_alf()
        self.set_alf()
        for atom in self:
            atom.calculate_features(self)
    
    @property
    def features(self):
        try:
            return self._features
        except AttributeError:
            self.calculate_features()
            self._features = [atom.features for atom in self]
            return self._features

    @property
    def nfeatures(self):
        return len(self.features[0])
    
    def __len__(self):
        return len(self._atoms)

    def __delitem__(self, i):
        del self._atoms[i]

    def __getitem__(self, i):
        return self._atoms[i]

    def __str__(self):
        return "\n".join([str(atom) for atom in self])

    def __repr__(self):
        return str(self)
    
    def __sub__(self, other):
        for i, atom in enumerate(self):
            for jatom in other:
                if jatom == atom:
                    del self[i]
        return self


class Atom:
    """ Deals with 1 Atom / 1 Coordinate Line, self._atoms containts items of class Atom,
    they are printed nicely because of the __str__ method of Atom class"""
    ang2bohr = 1.88971616463
    counter = it.count(1)
    def __init__(self, coordinate_line):
        self.atom_type = ""
        self.atom_number = next(Atom.counter)

        self.x = 0
        self.y = 0
        self.z = 0

        self.read_input(coordinate_line)

        self._bonds = []
        self.__level = it.count(0)

        self.x_axis = None
        self.xy_plane = None

        self.features = []

    def read_input(self, coordinate_line):
        find_atom = coordinate_line.split()
        self.atom_type = find_atom[0]
        coordinate_line = next(re.finditer(r"(\s*[+-]?\d+.\d+([Ee][+-]?\d+)?){3}", coordinate_line)).group()
        coordinate_line = re.finditer(r"[+-]?\d+.\d+([Ee][+-]?\d+)?", coordinate_line)
        self.x = float(next(coordinate_line).group())
        self.y = float(next(coordinate_line).group())
        self.z = float(next(coordinate_line).group())

    def dist(self, other):
        dist = 0
        for icoord, jcoord in zip(self.coordinates, other.coordinates):
            dist += (icoord - jcoord)**2
        return np.sqrt(dist)

    def xdiff(self, other):
        return other.x - self.x
    
    def ydiff(self, other):
        return other.y - self.y

    def zdiff(self, other):
        return other.z - self.z

    def angle(self, atom1, atom2):
        temp = self.xdiff(atom1) * self.xdiff(atom2) + \
                self.ydiff(atom1) * self.ydiff(atom2) + \
                self.zdiff(atom1) * self.zdiff(atom2)
        return math.acos((temp / (self.dist(atom1) * self.dist(atom2))))

    def set_bond(self, jatom):
        if not jatom in self._bonds:
            self._bonds.append(jatom)

    @property
    def priority(self):
        level = next(self.__level)
        return self.get_priorty(level)

    def get_priorty(self, level):
        atoms = Atoms(self)
        for _ in range(level):
            atoms_to_add = []
            for atom in atoms:
                for bonded_atom in atom.bonds:
                    if not bonded_atom in atoms:
                        atoms_to_add.append(bonded_atom)
            atoms.add(atoms_to_add)
        return atoms.priority

    def reset_level(self):
        self.__level = it.count(0)

    def add_alf_atom(self, atom):
        if self.x_axis is None:
            self.x_axis = atom
        else:
            self.xy_plane = atom

    def C_1k(self):
        return [self.xdiff(self.x_axis) / self.dist(self.x_axis),
                self.ydiff(self.x_axis) / self.dist(self.x_axis),
                self.zdiff(self.x_axis) / self.dist(self.x_axis)]

    def C_2k(self):
        xdiff1 = self.xdiff(self.x_axis)
        ydiff1 = self.ydiff(self.x_axis)
        zdiff1 = self.zdiff(self.x_axis)

        xdiff2 = self.xdiff(self.xy_plane)
        ydiff2 = self.ydiff(self.xy_plane)
        zdiff2 = self.zdiff(self.xy_plane)

        sigma_fflux = -(xdiff1 * xdiff2 + ydiff1 * ydiff2 + zdiff1 * zdiff2) / (
                xdiff1 * xdiff1 + ydiff1 * ydiff1 + zdiff1 * zdiff1)

        y_vec1 = sigma_fflux * xdiff1 + xdiff2
        y_vec2 = sigma_fflux * ydiff1 + ydiff2
        y_vec3 = sigma_fflux * zdiff1 + zdiff2

        yy = math.sqrt(y_vec1 * y_vec1 + y_vec2 * y_vec2 + y_vec3 * y_vec3)

        y_vec1 /= yy
        y_vec2 /= yy
        y_vec3 /= yy

        return [y_vec1, y_vec2, y_vec3]

    def C_3k(self, C_1k, C_2k):
        xx3 = C_1k[1]*C_2k[2] - C_1k[2]*C_2k[1]
        yy3 = C_1k[2]*C_2k[0] - C_1k[0]*C_2k[2]
        zz3 = C_1k[0]*C_2k[1] - C_1k[1]*C_2k[0]

        return [xx3, yy3, zz3]

    def calculate_features(self, atoms, unit="bohr"):
        ang2bohr = Atom.ang2bohr
        if "ang" in unit.lower():
            ang2bohr = 1.0

        x_bond = self.dist(self.x_axis)
        xy_bond = self.dist(self.xy_plane)
        angle = self.angle(self.x_axis, self.xy_plane)

        self.features += [x_bond * ang2bohr]
        self.features += [xy_bond * ang2bohr]
        self.features += [angle]

        for jatom in atoms:
            if jatom.num in self.alf_nums:
                continue
            self.features += [self.dist(jatom) * ang2bohr]

            C_1k = self.C_1k()
            C_2k = self.C_2k()
            C_3k = self.C_3k(C_1k, C_2k)
            
            xdiff = self.xdiff(jatom)
            ydiff = self.ydiff(jatom)
            zdiff = self.zdiff(jatom)

            zeta1 = C_1k[0] * xdiff + C_1k[1] * ydiff + C_1k[2] * zdiff
            zeta2 = C_2k[0] * xdiff + C_2k[1] * ydiff + C_2k[2] * zdiff
            zeta3 = C_3k[0] * xdiff + C_3k[1] * ydiff + C_3k[2] * zdiff

            # Calculate Theta
            self.features += [math.acos(zeta3 / self.dist(jatom))]

            # Calculate Phi
            self.features += [math.atan2(zeta2, zeta1)]

    def to_angstroms(self):
        self.x /= Atom.ang2bohr
        self.y /= Atom.ang2bohr
        self.z /= Atom.ang2bohr
    
    def to_bohr(self):
        self.x *= Atom.ang2bohr
        self.y *= Atom.ang2bohr
        self.z *= Atom.ang2bohr

    @property
    def bonds(self):
        return Atoms(self._bonds)

    @property
    def mass(self):
        global type2mass
        return type2mass[self.atom_type]
    
    @property
    def radius(self):
        global type2rad
        return type2rad[self.atom_type]

    @property
    def coordinates(self):
        return [self.x, self.y, self.z]
    
    @property
    def atom_num(self):
        return f"{self.atom_type}{self.atom_number}"
    
    @property
    def name(self):
        return self.atom_num

    @property
    def atom(self):
        return f"{self.atom_type}"

    @property
    def atom_num_coordinates(self):
        return [self.atom_num] + self.coordinates
    
    @property
    def atom_coordinates(self):
        return [self.atom] + self.coordinates
    
    @property
    def coordinates_string(self):
        width = str(16)
        precision = str(8)
        return f"{self.x:{width}.{precision}f}{self.y:{width}.{precision}f}{self.z:{width}.{precision}f}"

    @property
    def num(self):
        return self.atom_number

    @property
    def type(self):
        return self.atom_type

    @property
    def alf(self):
        alf = [self]
        if not self.x_axis is None: alf.append(self.x_axis)
        if not self.xy_plane is None: alf.append(self.xy_plane)
        return alf

    @property
    def alf_nums(self):
        return [atom.num for atom in self.alf]

    def __str__(self):
        return f"{self.atom_type:<3s}{self.coordinates_string}"

    def __repr__(self):
        return str(self)
    
    def __eq__(self, other):
        return self.atom_num == other.atom_num

    def __hash__(self):
        return hash(str(self.num) + str(self.coordinates_string))



def features_and_atom_names(xyz_file):
    """ Returns features as 3D array, [atom][point][feature]
    Example: 10 points water xyz file would have shape (3, 10, 3) where 3 is the number of atoms,
    10 is the number of points, and 3 is the number of features"""

    trajectory = Trajectory(xyz_file, read=True)
    features = trajectory.features # features[point][atom][feature]
    features = np.swapaxes(features, 0, 1) # features[atom][point][feature]

    atom_names = trajectory.atom_names
    numbered_priorities = trajectory.priorities
    # map the numbered priorities to actual names of atoms
    atom_name_priorities = [list(map(lambda i: atom_names[i-1], prio)) for prio in numbered_priorities]
    for atom_name_prio in atom_name_priorities:
        del atom_name_prio[0] # removing central atom which is not going to be present in the 3D xyz array used in plotting
        # only non-central atoms are in this array in order x_axis atom, xyz_plane atom, followed by atoms that were in spherical polar coords
    
    # numbered_priorities =  [list(map(lambda i: i-1, prio)) for prio in numbered_priorities]
    # for numbered_prio in numbered_priorities:
    #     del numbered_prio[0]

    return features, atom_names, atom_name_priorities

########################################################################
#                     CREATING 4D ARRAY TO PLOT
# 
########################################################################

# could not install matplotlib due to package conflicts, this list consists of colors from matplitlib
# colors = ['#F0F8FF', '#FAEBD7', '#00FFFF', '#7FFFD4', '#F0FFFF', '#F5F5DC',
#     '#FFE4C4', '#000000', '#FFEBCD', '#0000FF', '#8A2BE2', '#A52A2A', '#DEB887',
#     '#5F9EA0', '#7FFF00', '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C', '#00FFFF',
#     '#00008B', '#008B8B', '#B8860B', '#A9A9A9', '#006400', '#A9A9A9', '#BDB76B', '#8B008B',
#     '#556B2F', '#FF8C00', '#9932CC', '#8B0000', '#E9967A', '#8FBC8F', '#483D8B', '#2F4F4F',
#     '#2F4F4F', '#00CED1', '#9400D3', '#FF1493', '#00BFFF', '#696969', '#696969', '#1E90FF', '#B22222',
#     '#FFFAF0', '#228B22', '#FF00FF', '#DCDCDC', '#F8F8FF', '#FFD700', '#DAA520', '#808080', '#008000', '#ADFF2F',
#     '#808080', '#F0FFF0', '#FF69B4', '#CD5C5C', '#4B0082', '#FFFFF0', '#F0E68C', '#E6E6FA', '#FFF0F5', '#7CFC00',
#     '#FFFACD', '#ADD8E6', '#F08080', '#E0FFFF', '#FAFAD2', '#D3D3D3', '#90EE90', '#D3D3D3', '#FFB6C1', '#FFA07A',
#     '#20B2AA', '#87CEFA', '#778899', '#778899', '#B0C4DE', '#FFFFE0', '#00FF00', '#32CD32', '#FAF0E6', '#FF00FF',
#     '#800000', '#66CDAA', '#0000CD', '#BA55D3', '#9370DB', '#3CB371', '#7B68EE', '#00FA9A', '#48D1CC', '#C71585',
#     '#191970', '#F5FFFA', '#FFE4E1', '#FFE4B5', '#FFDEAD', '#000080', '#FDF5E6', '#808000', '#6B8E23', '#FFA500',
#     '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98', '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F', '#FFC0CB',
#     '#DDA0DD', '#B0E0E6', '#800080', '#663399', '#FF0000', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#F4A460',
#     '#2E8B57', '#FFF5EE', '#A0522D', '#C0C0C0', '#87CEEB', '#6A5ACD', '#708090', '#708090', '#FFFAFA', '#00FF7F',
#     '#4682B4', '#D2B48C', '#008080', '#D8BFD8', '#FF6347', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFFFFF', '#F5F5F5',
#     '#FFFF00', '#9ACD32']

colors = ["red", "green", "blue", "orange", "purple", "pink"]

class XYZArrays:
    """ Class for converting to Cartesian space. 
    Creates a 3D array for each atom on which the ALF is centered (1 4D array total). Each 2D array in the 3D array
    consists of N_pointsx3 (because each point has x,y,z coords in 3D space) matrices, where each matrix contains
    the xyz coordinates of every atom that is NOT the atom on which the ALF is centered."""

    def __init__(self, all_atom_features, atom_names):

        self.all_atom_features = all_atom_features
        self.atom_names = atom_names
        self.n_atoms, self.n_points, self.n_features = all_atom_features.shape
        self.all_atom_4d_array = self.stack_one_atom_xyz_3D_arrays()

    def stack_one_atom_xyz_3D_arrays(self):
        """Iterates over all the atoms in the molecule. Every atom can be used as center for ALF. 
        Stacks together 3D array for xy plane atoms, as well as 3D array that defines rest of atoms 
        in xyz coordinates.

        Returns all_atom_4D_array , a 4D array of shape (N_atoms, N_atoms-1, N_points, 3)
        If an atom is the ALF center, all other atoms have to be represented as xyz coordinates (thus the N_atoms-1)
        Every point for all other atoms has x, y, and z coordinates (so this is there the 3 comes in)

        Example: Methanol has 6 atoms, when an atom is used as the ALF center only 5 atoms remain to be expressed
        as xyz coordinates. These atoms have n number of points (number of timesteps in trajectory .xyz file), with
        each point having x,y,z coordinates."""

        all_other_atom_3D_arrays = []

        for one_atom_features in all_atom_features:

            # xy_atom_3d_array, and polar_atoms_3d_array are both 3D atom arrays
            # (xy_atom_matrix is 2xN_pointsx3, and polar atom matrix
            # is of shape N_remaining_atomsxN_pointsx3 where N remaining atoms is N_atoms-3)

            xy_atom_3d_array = self.get_xy_plane_atom_3d_array(one_atom_features)
            polar_atoms_3d_array = self.get_polar_atom_3d_array(one_atom_features)

            # so now we stack these matrices into one 3D array that is the xyz coordinates 
            # for all atoms OTHER than the atom on which the ALF is centered,
            # shape ((2+N_remaining_atoms),N_points,3)
            one_atom_total_array = np.concatenate((xy_atom_3d_array, polar_atoms_3d_array), axis=0)
            all_other_atom_3D_arrays.append(one_atom_total_array)

        # finally we can stack these 3D arrays into one 4D array, which will contain all the info needed
        # for plotting every atom as the ALF center. This 4D array will be stored and then if can be used
        # to quickly remove/add atoms in the visualization
        all_atom_4d_array = np.stack([i for i in all_other_atom_3D_arrays])

        return all_atom_4d_array

    def get_xy_plane_atom_3d_array(self, one_atom):

        """ Input: Takes in one atom feature matrix.
        Takes first three features for one atom and gives xyz coordinates of the two atoms used to define the x-axis, and the
        xy plane respectively
        """

        # gives an n_pointx1 vector of bond 1 lengths, need to convert to matrix with xyz coordinates (n_pointsx3)
        # y and z dimension are always 0s 
        tmp_bond1 = one_atom[:,[0]]
        z = np.zeros((tmp_bond1.shape[0], 2), dtype=tmp_bond1.dtype)
        bond1 = np.concatenate((tmp_bond1,z), axis=1)

        tmp_bond2 = one_atom[:,[1]] # bond 2 length from features (n_pointsx1)
        angle12 = one_atom[:,[2]] # angle between bond1 and bond2 (n_pontsx1)
        x_bond2 = np.multiply(np.cos(angle12), tmp_bond2)
        y_bond2 = np.multiply(np.sin(angle12), tmp_bond2)
        # z direction is always 0
        z_bond2 = np.zeros((tmp_bond2.shape[0], 1), dtype=tmp_bond2.dtype)
        bond2 = np.concatenate((x_bond2,y_bond2,z_bond2), axis=1)

        bond1_bond2_matrix = np.stack((bond1, bond2), axis=0)

        return bond1_bond2_matrix

    def get_polar_atom_3d_array(self, one_atom):

        """ Input: Takes in one atom feature matrix. 
        Every three features (after the first three features that define the xy plane) have a radius, theta, and
        phi component that defines where an atom is in xyz coordinates
        This method will return a 3D array of shape n_remaining_atoms, n_points, 3, wher 3 is because x,y,z coordinate
        """

        # firist three features account for 3 atoms (central atom, atom that defines x axis, and atom that defines 
        # xy plane. Therefore do not have to iterate over them)
        n_remaining_atoms = self.n_atoms - 3
        i,j,k = 3, 4, 5
        xyz_atoms_list = []

        for _ in range(n_remaining_atoms):

            r_data = one_atom[:,[i]] # vector of r distances (n_pointsx1)
            theta_data = one_atom[:,[j]] # vector of thetas (n_pointsx1)
            phi_data = one_atom[:,[k]] # vector of phis (n_pointsx1)
            xx = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.cos(phi_data))
            yy = np.multiply(np.multiply(r_data,np.sin(theta_data)),np.sin(phi_data))
            zz = np.multiply(r_data, np.cos(theta_data))
            xyz_atom = np.concatenate((xx,yy,zz), axis=1) # an N_pointsx3 matrix (storing xyz info for one atom)
            # polar_atoms_xyz = np.dstack((polar_atoms_xyz, xyz_atom), axis=2)
            xyz_atoms_list.append(xyz_atom)

            i += 3
            j += 3
            k += 3

        polar_atoms_xyz_matrix = np.stack([i for i in xyz_atoms_list])

        return polar_atoms_xyz_matrix


#########################################################################
#                       PYVISTA/ QT PLOTTING TOOL
#########################################################################
Ui_MainWindow, Ui_BaseClass = uic.loadUiType("testing.ui")

class VisualizationWindowDecorators:

    @staticmethod
    def clear_plot_add_grid(original_method):

        def wrapper(instance_reference):

            instance_reference.plotter.clear()
            original_method(instance_reference)
            instance_reference.plotter.show_grid()

        return wrapper

class VisualizationWindow(Ui_BaseClass):
    """ handles GUI and connects user commands with what to plot on pyvista plot"""

    def __init__(self, all_atom_dict, atom_names, atom_colors):

        super().__init__()

        self.all_atom_dict = all_atom_dict
        self.atom_names = atom_names # list of atom names
        self.atom_colors = atom_colors #dict of atom:color

        # used to initialize UI to plot first central alf atom (based on index, ex. C1, O1, etc.)
        self.current_central_atom_name = atom_names[0]
        self.current_central_atom_color = self.atom_colors[self.current_central_atom_name]
        self.center = np.array([0, 0, 0])

        # keeps total noncentral data that can be plotted
        self.all_noncentral_data = all_atom_dict[self.current_central_atom_name]
        self.current_noncentral_data = copy(self.all_noncentral_data)
        self.current_datablock = pv.MultiBlock(self.all_noncentral_data)

        # used in initializing values for slider, atom selecter, and atom color parts
        self.slider_position = 0
        self.current_noncentral_atom_names = [name for name in self.atom_names if name != self.current_central_atom_name]
        self.checkboxes = []
        self.current_checked_atoms = []

        self.color_buttons = []
        self.color_button_labels = []

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self._start_alf_vis_ui()

    def _start_alf_vis_ui(self):
        """ Initializes pyvista plot and user ui, with first atom ALF displayed
        Methods starting with _ only called here"""

        # initialize ui values and plotter
        self._start_combo_central_atom_names()
        self._start_points_to_plot_slider()
        self._start_pyvista_plotter()
        self.update_central_atom_and_plot()
    
    def _start_combo_central_atom_names(self):
        """ method initializing atom names combo box from list of atom names"""

        self.ui.atom_names_combo.addItems(self.atom_names)
        self.ui.atom_names_combo.currentIndexChanged.connect(self.update_central_atom_and_plot)

    def _start_points_to_plot_slider(self):

        self.ui.points_slider.setMinimum(1)
        self.ui.points_slider.setMaximum(100)
        self.ui.points_slider.setTickInterval(2)
        self.ui.points_slider.valueChanged.connect(self.update_atom_data_and_plot)

    def _start_pyvista_plotter(self):
        """ method to initialize pyvista plot"""

        self.plotter = QtInteractor(self.ui.pyvista_frame)
        self.ui.horizontalLayout_3.addWidget(self.plotter.interactor)
        self.plotter.show_grid()

    @VisualizationWindowDecorators.clear_plot_add_grid
    def update_central_atom_and_plot(self):
        """ main method if a different central atom is chosen
        method runs ONLY when the central atom is changed in the GUI"""
        """ Updates central atom (always at 0,0,0 but can update color if different atom) as 
        well as updates non central atom data"""

        self.update_central_atom_data()
        self.points_slider_reset()
        self.update_noncentral_atoms_data()
        self.update_checkboxes_widget()
        self.update_atom_color_box_buttons()
        self.plot_updated_data()

    @VisualizationWindowDecorators.clear_plot_add_grid
    def update_atom_data_and_plot(self):
        """main method for updating data/colors for the current central atom being shown
        method runs ONLY when there are changes in atoms to be plotted (based on GUI checkboxes changes), or
        when there are changes in atom colors as selected in the GUI"""

        self.update_checked_atoms()
        self.update_noncentral_atoms_data()
        self.update_atom_color_box_buttons()
        self.plot_updated_data()

    def update_central_atom_data(self):
        """ method used to update the central ALF atom, depending on selected atom in combo box"""

        self.current_central_atom_name = self.ui.atom_names_combo.currentText()
        self.current_central_atom_color = self.atom_colors[self.current_central_atom_name]

    def update_noncentral_atoms_data(self):

        self.current_noncentral_atom_names = [name for name in self.atom_names if name != self.current_central_atom_name]
        self.all_noncentral_data = self.all_atom_dict[self.current_central_atom_name]
        
        self.slider_update_plotted_data()

        self.current_datablock = pv.MultiBlock(self.current_noncentral_data)

    def points_slider_reset(self):
        """ reset slider when changing to a different central atom"""

        self.ui.points_slider.setValue(1)

    def slider_update_plotted_data(self):
        """ update noncentral atom data to plot when slider is moved. This makes slices of the data before making 
        it a pyvista MultiBlock"""

        slicing = self.ui.points_slider.value()
        for atom in self.all_noncentral_data.keys():
            self.current_noncentral_data[atom] = self.all_noncentral_data[atom][0::slicing,:]

    def update_checkboxes_widget(self):

        if self.checkboxes != []:
            for check in self.checkboxes:
                self.ui.gridLayout.removeWidget(check)
                check.deleteLater()
                check = None
            self.checkboxes = []

        self.current_checked_atoms = []
        row = 0
        column = 0
        for atom in self.current_noncentral_atom_names:
            checkbox = QtWidgets.QCheckBox(f"{atom}")
            checkbox.setCheckState(QtCore.Qt.Checked)
            checkbox.stateChanged.connect(self.update_atom_data_and_plot)
            self.checkboxes.append(checkbox)
            self.current_checked_atoms.append(checkbox.text())
            self.ui.gridLayout.addWidget(checkbox, row, column)
            column += 1
            if column % 3 == 0:
                row += 1
                column = 0

    def update_checked_atoms(self):

        for checkbox in self.checkboxes:
            if checkbox.isChecked() == False and checkbox.text() in self.current_checked_atoms:
                self.current_checked_atoms.remove(checkbox.text())
            elif checkbox.isChecked() == True and checkbox.text() not in self.current_checked_atoms:
                self.current_checked_atoms.append(checkbox.text())

    def update_atom_color_box_buttons(self):
        """ updates atoms that are in the color box, depending on which central atom is chosen and also which
        checkboxes are ticked in the checkbox box."""

        # clear color buttons and labels if checkboxes are changed
        if self.color_buttons != []:
            for button, button_label in zip(self.color_buttons, self.color_button_labels):
                self.ui.gridLayout.removeWidget(button)
                self.ui.gridLayout.removeWidget(button_label)
                button.deleteLater()
                button_label.deleteLater()
                button = None
                button_label = None
            self.color_buttons = []
            self.color_button_labels = []

        row = 0
        column_button = 0
        column_label = 1

        # add color button for central atom as well
        push_button = QtWidgets.QPushButton(self.current_central_atom_name)
        push_button.setStyleSheet(f"background-color : {self.atom_colors[self.current_central_atom_name]}; color: {self.atom_colors[self.current_central_atom_name]};")
        push_button.clicked.connect(self.change_atom_color)
        push_button_label = QtWidgets.QLabel(self.current_central_atom_name)
        self.color_buttons.append(push_button)
        self.color_button_labels.append(push_button_label)
        self.ui.gridLayout_2.addWidget(push_button, row, column_button)
        self.ui.gridLayout_2.addWidget(push_button_label, row, column_label)

        column_button += 2
        column_label += 2

        for checkbox in self.checkboxes:
            if checkbox.isChecked() == True:

                push_button = QtWidgets.QPushButton(f"{checkbox.text()}")
                push_button.setStyleSheet(f"background-color : {self.atom_colors[checkbox.text()]}; color: {self.atom_colors[checkbox.text()]};")
                push_button.clicked.connect(self.change_atom_color)

                push_button_label = QtWidgets.QLabel(f"{checkbox.text()}")

                self.color_buttons.append(push_button)
                self.color_button_labels.append(push_button_label)

                self.ui.gridLayout_2.addWidget(push_button, row, column_button)
                self.ui.gridLayout_2.addWidget(push_button_label, row, column_label)

                column_button += 2
                column_label += 2

                if column_button % 3 == 0:
                    row += 1
                    column_button = 0
                    column_label = 1

    def change_atom_color(self):
        """ opens up color dialog and lets user select a new color for a particular atom """
        
        color = QtWidgets.QColorDialog.getColor() # this gives a QColor Object

        # find which button called this method and change its color as needed
        for push_button in self.color_buttons:
            if push_button == QtCore.QObject.sender(self):

                if push_button.text() == self.current_central_atom_name:

                    self.current_central_atom_color = color.name()
                    push_button.setStyleSheet("")
                    push_button.setStyleSheet(f"background-color : {color.name()}; color: {color.name()};")
                    self.atom_colors[f"{push_button.text()}"] = color.name()
                else:

                    push_button.setStyleSheet("")
                    push_button.setStyleSheet(f"background-color : {color.name()}; color: {color.name()};")
                    self.atom_colors[f"{push_button.text()}"] = color.name()

        self.plot_updated_data()

    def menu_lock(self):
        """ lock menu when choosing colors so that people cannot switch plotted atoms while color dialog is opened.
        prevents bugs"""
        pass

    def plot_updated_data(self):
        """ plots data after an update to central ALF atom"""

        center = pv.PolyData(self.center)
        self.plotter.add_mesh(center, color=self.current_central_atom_color, point_size=30, render_points_as_spheres=True)

        for block in self.current_datablock.keys():
            if block in self.current_checked_atoms:
                color = self.atom_colors.get(block)
                self.plotter.add_mesh(self.current_datablock[block], color=color, point_size=10, render_points_as_spheres=True)

if __name__ == "__main__":

    xyz_files = sorted(glob("*.xyz"))
    if len(xyz_files) == 1:
        xyz_file = xyz_files[0]
    else:
        print("Select xyz file to evaluate:")
        print ("")
        for i in range(len(xyz_files)):
            print (f"[{i+1}] --> ", xyz_files[i])
        xyz_file = xyz_files[int(input())-1]

    # all_atom_features are 3D array [atom][point][feature], shape is (n_atoms, n_points, n_features)
    all_atom_features, atom_names, atoms_names_priorities = features_and_atom_names(xyz_file)
    atom_colors = dict(zip(atom_names, cycle(colors))) # initialize atom colors

    system_as_xyz = XYZArrays(all_atom_features, atom_names) # gives 4D numpy array

    total_dict = {} # dictionary of dictionaries, ordered as {"C1":{"O3":xyz_array, "H2":xyz_array, "H4":xyz_array ....}, 
    # "H2": {"C1":xyz_array, "O3":xyz_array.......}........,"H6":{"O3":xyz_array, "C1":xyz_array}
    # the ordering is {central_atom1: {x_axis_atom:xyz_array, xy_plane_atom:xyz_array, polar_atom:xyz_array ..... },
    #  central_atom2:{x_axis_atom:xyz_array, xy_plane_atom:xyz_array, polar_atom:xyz_array, ..... }....}
    # this was done to keep track of colors (which cannot really be done using np arrays)

    xyz_dict = {} # this dict keeps inner dictionaries from the total_array such as {"O3":xyz_array, "H2":xyz_array, "H4":xyz_array ....}
    # it gets reset after every iteration of the loop to move onto the next atom center

    # iterate over central atoms, their respective 3D array of other atoms as xyz coords, as well as the priorities of these xyz atoms (i.e which is
    # x axis, which is xy plane etc.)
    for center_atom, center_atom_xyzs, atom_names_prio in zip(atom_names, system_as_xyz.all_atom_4d_array, atoms_names_priorities):
        # C1  #C1 non central atom 3D array # O3, H2, H4, H5, H6 # 2, 1, 3, 4, 5
        for idx, atom_name_prio in enumerate(atom_names_prio):
            #  0 O3, 1 H2, etc. used to get out the individual 2d arrays which are ordered as 
            # x axis, xy plane, and then all atoms in polar coords
            xyz_dict[atom_name_prio] = center_atom_xyzs[idx]
        
        total_dict[center_atom] = xyz_dict
        xyz_dict = {}

    app = QtWidgets.QApplication(sys.argv)
    main_window = VisualizationWindow(total_dict, atom_names, atom_colors)
    main_window.show()

    app.exec_()