#!/usr/bin/env python

import os
import pybel


class LPC:
    def __init__(self, lig_path, prt_path):
        self._lig_path = lig_path
        self._prt_path = prt_path
        self._lig_format = os.path.splitext(self._lig_path)[-1][1:]
        self._prt_format = os.path.splitext(self._prt_path)[-1][1:]

    @staticmethod
    def merge(prt, lig):
        prt_pdb_lines = filter(lambda s: 'ATOM' in s,
                               prt.write('pdb').splitlines(True))
        lig_pdb_lines = filter(lambda s: ('ATOM' in s) or ('HETATM' in s),
                               lig.write('pdb').splitlines(True))
        to_write = []
        to_write.append("MODEL 1\n")
        to_write.extend(prt_pdb_lines)
        to_write.append("TER\n")
        to_write.extend(lig_pdb_lines)
        to_write.append("END\n")
        return "".join(to_write)

    def runLPC(self):
        prt = pybel.readfile(self._prt_format, self._prt_path).next()
        lig = pybel.readfile(self._lig_format, self._lig_path).next()
        self.merge(prt, lig)


def test():
    lpc = LPC("/work/jaydy/dat/edud/crystal-ligands-sdf/1b9v.sdf",
              "/work/jaydy/dat/edud/crystal-structure/1b9vA.pdb")
    lpc.runLPC()


if __name__ == '__main__':
    test()
