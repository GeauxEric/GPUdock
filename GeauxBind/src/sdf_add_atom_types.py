import pybel


authentic_sdf = "../data/1a07C1.sdf"


with open(authentic_sdf, 'r') as ifs:
    lines = ifs.readlines()
    for idx, line in enumerate(lines):
        if "OB_ATOM_TYPES" in line:
            types = lines[idx + 1].rstrip()


lig = pybel.readfile("sdf", authentic_sdf).next()

types_from_ob = []
for atom in lig:
    print atom.type
    types_from_ob.append(atom.type)

assert ' '.join(types_from_ob) == types
