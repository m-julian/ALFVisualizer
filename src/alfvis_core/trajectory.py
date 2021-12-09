import re
from pathlib import Path
from typing import List
import numpy as np
from alfvis_core.list_of_atoms import ListOfAtoms
from alfvis_core.atoms import Atoms
from alfvis_core.atom import Atom
from alfvis_core.file_state import FileState
import ast


class Trajectory(ListOfAtoms):
    """Handles .xyz files that have multiple timesteps, with each timestep giving the x y z coordinates of the
    atoms."""

    def __init__(self, path: Path):
        ListOfAtoms.__init__(self)
        self.state = FileState.Unread
        self.path = path

    def _read_file(self):
        with open(self.path, "r") as f:
            atoms = Atoms()
            for line in f:
                if not line.strip():
                    continue
                elif re.match(r"^\s*\d+$", line):
                    natoms = int(line)
                    continue
                # this is the comment line that can contain extra info
                elif re.match(r"^\s*?i\s*?=\s*?\d+\s*energy", line):
                    properties_error = line.split("=")[-1].strip()
                    atoms.properties_error = ast.literal_eval(properties_error)
                while len(atoms) < natoms:
                    line = next(f)
                    if re.match(
                        r"\s*?\w+(\s+[+-]?\d+.\d+([Ee]?[+-]?\d+)?){3}", line
                    ):
                        atom_type, x, y, z = line.split()
                        atoms.add(
                            Atom(atom_type, float(x), float(y), float(z))
                        )
                self.add(atoms)
                atoms = Atoms()

        self.state = FileState.Read

    @property
    def filetype(self) -> str:
        return ".xyz"

    def add(self, atoms):
        """Add a list of Atoms (corresponding to one timestep) to the end of the trajectory list"""
        if isinstance(atoms, Atoms):
            self.append(atoms)
        else:
            self.append(Atoms(atoms))

    def extend(self, atoms):
        """extend the current trajectory by another list of atoms (could be another trajectory list)"""
        if isinstance(atoms, Atoms):
            self.extend(atoms)
        else:
            self.extend(Atoms(atoms))

    def write(self, fname=None):
        """write a new .xyz file that contains the timestep i, as well as the coordinates of the atoms
        for that timestep."""
        if fname is None:
            fname = self.path
        with open(fname, "w") as f:
            for i, atoms in enumerate(self):
                f.write(f"    {len(atoms)}\ni = {i}\n")
                f.write(f"{atoms}\n")

    def rmsd(self, ref=None):
        if ref is None:
            ref = self[0]
        elif isinstance(ref, int):
            ref = self[ref]

        return [ref.rmsd(point) for point in self]

    @property
    def properties_error(self) -> List:
        """returns a list of energies as read in from the .xyz file comment line.
        This is used to plot colormaps of the whole trajectory to see any points which produce poor results."""
        if hasattr(self[0], 'properties_error'):
            return [timesteps.properties_error for timesteps in self]
        # if no energy/errors have been read in this needs to return None because it is used in the GUI class __int__
        else:
            return None

    @property
    def features(self) -> np.ndarray:
        """
        Returns:
            :type: `np.ndarray`
            A 3D array of features for every atom in every timestep. Shape `n_timesteps` x `n_atoms` x `n_features`)
            If the trajectory instance is indexed by str, the array has shape `n_timesteps` x `n_features`.
            If the trajectory instance is indexed by str, the array has shape `n_atoms` x `n_features`.
            If the trajectory instance is indexed by slice, the array has shape `slice`, `n_atoms` x `n_features`.
        """
        return np.array([timestep.features for timestep in self])

    @property
    def coordinates(self) -> np.ndarray:
        """
        Returns:
            :type: `np.ndarray`
            the xyz coordinates of all atoms for all timesteps. Shape `n_timesteps` x `n_atoms` x `3`
        """
        return np.array([timestep.coordinates for timestep in self])

    @property
    def alf(self):
        """ Returns the Atomic Local Frame (ALF) of the first Atoms instance in the trajectory. We expect the ALF to remain the same for all atoms
        throughout the trajectory. Atoms are indexed as in Python lists (start at 0)"""
        return self[0].alf

    @property
    def alf_index(self):
        """ Returns the Atomic Local Frame (ALF) of the first Atoms instance in the trajectory. We expect the ALF to remain the same for all atoms
        throughout the trajectory. Atoms are indexes as in their names (start at 1)"""
        return self[0].alf_index

    @property
    def priorities(self):
        alf_list = self.alf.tolist()

        for alf_atom in alf_list:
            for n_atom in range(0, len(self[0])):
                if n_atom in alf_atom:
                    continue
                else:
                    alf_atom.append(n_atom)

        return alf_list

    def __getitem__(self, item):
        # TODO: Implement isinstance(item, slice) if needed
        """ Used to index a Trajectory instance by a str (eg. trajectory['C1']) or by integer (eg. trajectory[2]),
        remember that indeces in Python start at 0, so trajectory[2] is the 3rd timestep.
        You can use something like (np.array([traj[i].features for i in range(2)]).shape) to features of a slice of
        a trajectory as slice is not implemented in __getitem__"""
        if self.state is not FileState.Read:
            self._read_file()

        return super().__getitem__(item)

    def __iter__(self):
        """ Used to iterate over timesteps (Atoms instances) in places such as for loops"""
        if self.state is not FileState.Read:
            self._read_file()
        return super().__iter__()

    def __len__(self):
        """ Returns the number of timesteps in the Trajectory instance"""
        if self.state is not FileState.Read:
            self._read_file()
        return super().__len__()
