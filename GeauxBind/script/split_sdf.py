#!/usr/bin/env python

import os


def split(sdf_path, dirname=""):
    """split the molecules in one sdf into individual files
    Keyword Arguments:
    sdf_path -- file path
    dirname  -- output directory, the same directory as the input file by default
    """
    section = []
    if dirname == "":
        dirname = os.path.dirname(sdf_path)
    ofn = ""
    with open(sdf_path, 'r') as ifs:
        lines = ifs.readlines()
        for idx, line in enumerate(lines):
            section.append(line)
            if '>  <MOLID>' in line:
                molecule_name = lines[idx + 1].rstrip()
                ofn = os.path.join(dirname, molecule_name + '.sdf')
            if '$$$$' in line:
                with open(ofn, 'w') as ofs:
                    ofs.write("".join(section))
                section = []


def main():
    split("/work/jaydy/tmp/1b9vA_sdf/1b9vA_3.sdf")
    pass


if __name__ == '__main__':
    main()
