#
# TAKES .XYZ TRAJECTORY FILE AND CALCULATES FEATURES OF EACH ATOM IN EACH POINT 
#

import re
import numpy as np
import itertools as it
import math
from functools import lru_cache
import xlsxwriter
import string
import pandas as pd
from glob import glob

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
    def __init__(self, fname, read=False):
        self.fname = fname
        self._trajectory = []

        if self.fname and read:
            self.read()
    
    def read(self):
        with open(self.fname, "r") as f:
            atoms = Atoms()
            for line in f: # reads through lines in .xyz file
                if not line.strip():
                    continue # ignore blank lines
                elif re.match(r"^\s+\d+$", line): # finds number of atoms in xyz file
                    natoms = int(line)
                    while len(atoms) < natoms:
                        line = next(f)
                        if re.match(r"\s*\w+(\s+[+-]?\d+.\d+([Ee]?[+-]?\d+)?){3}", line): # if re finds atom type and x, y,z coordinates it will add to list
                            atoms.add(line)
                    self.add(atoms)
                    atoms = Atoms()


    def append(self, atoms):
        self._trajectory.append(atoms)
    
    def add(self, atoms):
        self.append(atoms)

    @property
    def features(self):
        return np.array([point.features for point in self])

    @property
    def atom_names(self):
        return [atom.name for atom in self[0]]

    @property
    def nfeatures(self):
        return self[0].nfeatures

    def __len__(self):
        return len(self._trajectory)
    
    def __getitem__(self, i):
        return self._trajectory[i]

class Atom:
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
    def priority(self):
        level = next(self.__level)
        return self.get_priorty(level)

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


class Atoms:
    ALF = []

    def __init__(self, atoms=None):
        Atom.counter = it.count(1)
        self._atoms = []
        self._connectivity = None

        if not atoms is None:
            self.add(atoms)

    def add(self, atom):
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

    def __len__(self):
        return len(self._atoms)

    def __delitem__(self, i):
        del self._atoms[i]

    def __getitem__(self, i):
        return self._atoms[i]
    
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

    def __str__(self):
        return "\n".join([str(atom) for atom in self])

    def __repr__(self):
        return str(self)
    
    def __sub__(self, other):
        # other = sorted(Atoms(other), key=lambda x: x.num, reverse=False)
        for i, atom in enumerate(self):
            for jatom in other:
                if jatom == atom:
                    del self[i]
        return self
    
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
        self.set_alf()

    def set_alf(self):
        # self._atoms = sorted(self._atoms, key=lambda x: x.num)
        for atom, atom_alf in zip(self, Atoms.ALF):
            atom.x_axis = self[atom_alf[1]-1]
            atom.xy_plane = self[atom_alf[2]-1]

    def calculate_features(self):
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



def get_headings(nfeatures):
    "writes headings for columns"

    headings = [
        # "Training Set Point", 
        # "Atom Label", 
        "bond1", 
        "bond2", 
        "angle"
        ]

    nfeatures = len(trajectory.features[0][0])-3 # Removes bond1, bond 2, angle
    reduced_nfeatures = int(nfeatures/3) # each feature has r, theta, phi component

    for nfeature in range(1,(reduced_nfeatures+1)): #range goes to i-1
        feature_distance = f"r{nfeature+2}"  # starts at r3
        feature_theta = f"theta{nfeature+2}"   #starts at theta3
        feature_phi = f"phi{nfeature+2}"   #starts at phi3
        headings.append(feature_distance)
        headings.append(feature_theta)
        headings.append(feature_phi)

    return headings

def features_and_atom_names(xyz_file):
    """ Returns features as 3D array, [atom][point][feature]
    Example: 10 points water xyz file would have shape (3, 10, 3) where 3 is the number of atoms,
    10 is the number of points, and 3 is the number of features"""

    trajectory = Trajectory(xyz_file, read=True)
    features = trajectory.features # features[point][atom][feature]
    features = np.swapaxes(features, 0, 1) # features[atom][point][feature]
    atom_names = trajectory.atom_names

    return features, atom_names


