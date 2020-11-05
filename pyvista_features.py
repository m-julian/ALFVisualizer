import calculate_alf
from glob import glob
import pyvista as pv
from itertools import cycle
import numpy as np

# could not install matplotlib due to package conflicts, this list consists of colors from matplitlib
colors = ['#F0F8FF', '#FAEBD7', '#00FFFF', '#7FFFD4', '#F0FFFF', '#F5F5DC',
    '#FFE4C4', '#000000', '#FFEBCD', '#0000FF', '#8A2BE2', '#A52A2A', '#DEB887',
    '#5F9EA0', '#7FFF00', '#D2691E', '#FF7F50', '#6495ED', '#FFF8DC', '#DC143C', '#00FFFF',
    '#00008B', '#008B8B', '#B8860B', '#A9A9A9', '#006400', '#A9A9A9', '#BDB76B', '#8B008B',
    '#556B2F', '#FF8C00', '#9932CC', '#8B0000', '#E9967A', '#8FBC8F', '#483D8B', '#2F4F4F',
    '#2F4F4F', '#00CED1', '#9400D3', '#FF1493', '#00BFFF', '#696969', '#696969', '#1E90FF', '#B22222',
    '#FFFAF0', '#228B22', '#FF00FF', '#DCDCDC', '#F8F8FF', '#FFD700', '#DAA520', '#808080', '#008000', '#ADFF2F',
    '#808080', '#F0FFF0', '#FF69B4', '#CD5C5C', '#4B0082', '#FFFFF0', '#F0E68C', '#E6E6FA', '#FFF0F5', '#7CFC00',
    '#FFFACD', '#ADD8E6', '#F08080', '#E0FFFF', '#FAFAD2', '#D3D3D3', '#90EE90', '#D3D3D3', '#FFB6C1', '#FFA07A',
    '#20B2AA', '#87CEFA', '#778899', '#778899', '#B0C4DE', '#FFFFE0', '#00FF00', '#32CD32', '#FAF0E6', '#FF00FF',
    '#800000', '#66CDAA', '#0000CD', '#BA55D3', '#9370DB', '#3CB371', '#7B68EE', '#00FA9A', '#48D1CC', '#C71585',
    '#191970', '#F5FFFA', '#FFE4E1', '#FFE4B5', '#FFDEAD', '#000080', '#FDF5E6', '#808000', '#6B8E23', '#FFA500',
    '#FF4500', '#DA70D6', '#EEE8AA', '#98FB98', '#AFEEEE', '#DB7093', '#FFEFD5', '#FFDAB9', '#CD853F', '#FFC0CB',
    '#DDA0DD', '#B0E0E6', '#800080', '#663399', '#FF0000', '#BC8F8F', '#4169E1', '#8B4513', '#FA8072', '#F4A460',
    '#2E8B57', '#FFF5EE', '#A0522D', '#C0C0C0', '#87CEEB', '#6A5ACD', '#708090', '#708090', '#FFFAFA', '#00FF7F',
    '#4682B4', '#D2B48C', '#008080', '#D8BFD8', '#FF6347', '#40E0D0', '#EE82EE', '#F5DEB3', '#FFFFFF', '#F5F5F5',
    '#FFFF00', '#9ACD32']

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
all_atom_features, atom_names = calculate_alf.features_and_atom_names(xyz_file)
atom_colors = dict(zip(atom_names, cycle(colors))) # initialize atom colors

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
        i,j,k = 3, 3, 3
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


class GuiTool:
    """ Class handling dearpygui user interface"""

    def __init__(self, atom_names, current_central_atom):
        """ atom_names: list of atom names in molecule"""

        self.atom_names = atom_names
        self.n_atoms = len(atom_names)
        self.current_central_atom = atom_names[0] # inialize to first atom in molecule
        self.atoms_to_show = atom_names.remove(self.current_central_atom)

    def switch_central_atom(self, new_central_atom):

        self.current_central_atom = new_central_atom

    def remove_unwanted_atoms(self, atoms_to_remove):

        pass

class PlottingTool:
    """ Class for pyvista visualizaton of ALF"""

    def __init__(self, center_atom, atoms_to_plot, all_atom_4d_array):
        """ center_atom: atom on which ALF is centered
        atoms_to_plot: atoms to plot in pyvista visualization
        all_atom_4d_array: the 4d array created in XYZArrays class"""

        self.center_atom = center_atom
        self.atoms_to_plot = atoms_to_plot

if __name__ == "__main__":

    system_as_xyz = XYZArrays(all_atom_features, atom_names)

    print(system_as_xyz.all_atom_4d_array)
