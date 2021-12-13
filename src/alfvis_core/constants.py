from typing import Dict, List
import numpy as np

# colors

random_colors = ['#00FFFF', '#7FFFD4', '#FFE4C4', '#000000', '#DEB887',
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

# use these default colors for each atom type if you do not care about looking at specific atoms
default_atom_colors = {"O": "red", "H": "white", "Cl": "green", "N": "blue", "C": "grey", "S": "yellow", "P": "orange", "F": "purple"}


type2mass: Dict[str, float] = {
    "H": 1.007825,
    "He": 4.002603,
    "Li": 7.016005,
    "Be": 9.012182,
    "B": 11.009305,
    "C": 12.0,
    "N": 14.003074,
    "O": 15.994915,
    "F": 18.998403,
    "Ne": 19.99244,
    "Na": 22.989769,
    "Mg": 23.985042,
    "Al": 26.981539,
    "Si": 27.976927,
    "P": 30.973762,
    "S": 31.972071,
    "Cl": 34.968853,
    "Ar": 39.962383,
    "K": 38.963707,
    "Ca": 39.962591,
    "Sc": 44.955912,
    "Ti": 47.947946,
    "V": 50.94396,
    "Cr": 51.940508,
    "Mn": 54.938045,
    "Fe": 55.9349382,
    "Co": 58.933195,
    "Ni": 57.935343,
    "Cu": 62.929598,
    "Zn": 63.929142,
    "Ga": 68.925574,
    "Ge": 73.921178,
    "As": 74.921597,
    "Se": 79.916521,
    "Br": 78.918337,
    "Kr": 83.911507,
}

type2rad: Dict[str, float] = {
    "H": 0.37,
    "He": 0.32,
    "Li": 1.34,
    "Be": 0.9,
    "B": 0.82,
    "C": 0.77,
    "N": 0.74,
    "O": 0.73,
    "F": 0.71,
    "Ne": 0.69,
    "Na": 1.54,
    "Mg": 1.3,
    "Al": 1.18,
    "Si": 1.11,
    "P": 1.06,
    "S": 1.02,
    "Cl": 0.99,
    "Ar": 0.97,
    "K": 1.96,
    "Ca": 1.74,
    "Sc": 1.44,
    "Ti": 1.36,
    "V": 1.25,
    "Cr": 1.27,
    "Mn": 1.39,
    "Fe": 1.25,
    "Co": 1.26,
    "Ni": 1.21,
    "Cu": 1.38,
    "Zn": 1.31,
    "Ga": 1.26,
    "Ge": 1.22,
    "As": 1.19,
    "Se": 1.16,
    "Br": 1.14,
    "Kr": 1.1,
}

type2valence: Dict[str, int] = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "Ne": 8,
    "Na": 9,
    "Mg": 10,
    "Al": 3,
    "Si": 4,
    "P": 5,
    "S": 6,
    "Cl": 7,
    "Ar": 8,
    "K": 9,
    "Ca": 10,
    "Sc": 11,
    "Ti": 12,
    "V": 13,
    "Cr": 14,
    "Mn": 15,
    "Fe": 16,
    "Co": 17,
    "Ni": 18,
    "Cu": 11,
    "Zn": 12,
    "Ga": 13,
    "Ge": 4,
    "As": 5,
    "Se": 6,
    "Br": 7,
    "Kr": 8,
    "Sr": 10,
    "Y": 11,
    "Zr": 12,
    "Mo": 14,
    "Ru": 16,
    "Rh": 17,
    "Pd": 18,
    "Ag": 11,
    "In": 13,
    "Sb": 5,
    "Te": 6,
    "I": 7,
    "Ba": 10,
    "Ce": 12,
    "Gd": 18,
    "W": 14,
    "Au": 11,
    "Bi": 5,
}

dlpoly_weights: Dict[str, float] = {
    "H": 1.007975,
    "He": 4.002602,
    "Li": 6.9675,
    "Be": 9.0121831,
    "B": 10.8135,
    "C": 12.0106,
    "N": 14.006855,
    "O": 15.9994,
    "F": 18.99840316,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3055,
    "Al": 26.9815385,
    "Si": 28.085,
    "P": 30.973762,
    "S": 32.0675,
    "Cl": 35.4515,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955908,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938044,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.63,
    "As": 74.921595,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90584,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Mo": 95.95,
}

# ha_to_kj_mol = 2625.5
ha_to_kj_mol: float = (
    2625.4996394799  # Taken from https://en.wikipedia.org/wiki/Hartree
)
# The wikipedia article is converted from https://physics.nist.gov/cgi-bin/cuu/Value?hr

bohr2ang = 0.529177210903  # Converted from https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
ang2bohr = 1.0 / bohr2ang

rt3: float = np.sqrt(3)
rt5: float = np.sqrt(5)
rt6: float = np.sqrt(6)
rt10: float = np.sqrt(10)
rt15: float = np.sqrt(15)
rt35: float = np.sqrt(35)
rt70: float = np.sqrt(70)
rt1_24: float = np.sqrt(1 / 24)
rt_1_5: float = np.sqrt(1 / 5)
rt_1_10: float = np.sqrt(1 / 10)
rt_1_35: float = np.sqrt(1 / 35)
rt2_3: float = np.sqrt(2 / 3)
rt_2_35: float = np.sqrt(2 / 35)
rt3_3: float = np.sqrt(3) / 3
rt_3_3: float = np.sqrt(3) / 2
rt_3_5: float = np.sqrt(3 / 5)
rt3_8: float = np.sqrt(3 / 8)
rt5_8: float = np.sqrt(5 / 8)
rt5_12: float = np.sqrt(5 / 12)
rt_8_5: float = np.sqrt(8 / 5)
rt12_3: float = np.sqrt(12) / 3