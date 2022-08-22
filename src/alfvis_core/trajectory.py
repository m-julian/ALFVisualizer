from ichor.core.files import Trajectory
from ichor.core.atoms import Atom, Atoms
import re
from typing import List
import ast

class ALFVisTrajectory(Trajectory):

    def __init__(self, path, *args, **kwargs):

        super().__init__(path, *args, **kwargs)

    def _read_file(self):

        with open(self.path, "r") as f:

            # make empty Atoms instance in which to store one timestep
            atoms = Atoms()

            for line in f:

                # match the line containing the number of atoms in timestep
                if re.match(r"^\s*\d+$", line):
                    natoms = int(line)

                    # this is the comment line of xyz files. It can be empty or contain some useful information that can be stored.
                    line = next(f)

                    # if the comment line properties errors, we can store these
                    if re.match(r"^\s*?i\s*?=\s*?\d+\s*properties_error", line):
                        properties_error = line.split("=")[-1].strip()
                        atoms.properties = ast.literal_eval(properties_error)

                    # the next line after the comment line is where coordinates begin
                    for _ in range(natoms):
                        line = next(f)
                        if re.match(r"^\s*?\w+(\s+[+-]?\d+.\d+([Ee]?[+-]?\d+)?){3}", line):
                            atom_type, x, y, z = line.split()
                            atoms.add(Atom(atom_type, float(x), float(y), float(z)))

                    # add the Atoms instance to the Trajectory instance
                    self.add(atoms)
                    # make new Atoms instance where next timestep can be stored
                    atoms = Atoms()

    @property
    def properties(self) -> List:
        """returns a list of energies as read in from the .xyz file comment line.
        This is used to plot colormaps of the whole trajectory to see any points which produce poor results."""
        if hasattr(self[0], "properties_error"):
            return [timesteps.properties for timesteps in self]
        # if no energy/errors have been read in this needs to return None because it is used in the GUI class __int__
        else:
            return None
