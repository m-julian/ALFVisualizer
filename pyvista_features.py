import calculate_alf
from glob import glob

xyz_files = sorted(glob("*.xyz"))
if len(xyz_files) == 1:
    xyz_file = xyz_files[0]
else:
    print("Select xyz file to evaluate:")
    print ("")
    for i in range(len(xyz_files)):
        print (f"[{i+1}] --> ", xyz_files[i])
    xyz_file = xyz_files[int(input())-1]

features, atom_names = calculate_alf.features_and_atom_names(xyz_file)

print(atom_names)