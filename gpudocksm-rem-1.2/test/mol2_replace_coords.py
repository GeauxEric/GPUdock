################################################################################
ifn_mol = '1a07C1.mol2'
ifn_coords = 'LigCoord.txt'
ofn_mol = 'new_1a07C1.mol2'

mol_f = open(ifn_mol).readlines()
coord_f = open(ifn_coords).readlines()

begin, end = 0, 0
for i, line in enumerate(coord_f):
    if "@<BEGIN>ATOM" in line:
        begin = i
    if "@<END>ATOM" in line:
        end = i

lna = end - begin - 1
coords = coord_f[begin+1: begin+lna+1]
coords = [line.rstrip('\n') for line in coords]
coords = [line.split() for line in coords]


mol_begin, mol_end = 0, 0
for i, line in enumerate(mol_f):
    if '@<TRIPOS>ATOM' in line:
        mol_begin = i
    if '@<TRIPOS>BOND' in line:
        mol_end = i
old_coord_lines = mol_f[mol_begin+1:mol_end]

def replaceXYZ(old_line, new_xyz):
    """
    replace the x y z coords in the old file with the new x y z
    old x y z supposed to be indexed at 2:5
    each coordinate takes 10 digit
    """
    old_xyz = old_line.split()[2:5]
    old_z = old_xyz[2]
    right_end = old_line.index(old_z) + len(old_z)
    
    x, y, z = new_xyz
    # each coordinate takes 10 digit
    new_sub_str = '{:>10}{:>10}{:>10}'.format(x, y, z)

    old_line_lst = list(old_line)
    old_line_lst[right_end-30:right_end] = new_sub_str
    
    new_line = ''.join(old_line_lst)

    return new_line

new_coord_lines = [replaceXYZ(old_line, coords[idx]) for idx, old_line in enumerate(old_coord_lines)]

mol_f[mol_begin+1:mol_end] = new_coord_lines


with open(ofn_mol, 'w') as f:
    f.writelines(mol_f)
