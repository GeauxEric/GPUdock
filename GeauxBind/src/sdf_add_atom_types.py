#!/usr/bin/env python

import pybel
from openbabel import OBTypeTable

typetable = OBTypeTable()
typetable.SetFromType('INT')
typetable.SetToType('SYB')


def addTypes(ifn, ofn):
    """add the SYB atom types to the molecules in a sdf file
    """
    largeSDfile = pybel.Outputfile("sdf", ofn)

    for lig in pybel.readfile("sdf", ifn):
        types_from_ob = []
        for atom in lig:
            types_from_ob.append(typetable.Translate(atom.type))

        lig.data['NEW_OB_ATOM_TYPES'] = ' '.join(types_from_ob)
        largeSDfile.write(lig)

    largeSDfile.close()


def test():
    ifn = "../data/ZINC_3.sdf"
    ofn = ifn + ".add.sdf"
    addTypes(ifn, ofn)


if __name__ == '__main__':
    import sys
    ifn = sys.argv[1]
    ofn = sys.argv[2]
    addTypes(ifn, ofn)
