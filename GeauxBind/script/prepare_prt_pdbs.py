#!/usr/bin/env python

import luigi
import os
import shutil
from edud4geauxdock import Extract
from lpc import LPC
from prt_ens import genProteinEnsemble


class PrepareProteinEns(Extract):
    def output(self):
        path = os.path.join(self.dirname, self.ligand_code + '.pdb')
        return luigi.LocalTarget(path)

    @property
    def dirname(self):
        path = os.path.join(self.work_dir, 'protein_pdb')
        try:
            os.makedirs(path)
        except:
            pass
        return path

    @property
    def native_lig_sdf(self):
        return os.path.join(self.crystal_ligands_sdf, self.prt_code + '.sdf')

    @property
    def native_prt_pdb(self):
        return os.path.join(self.crystal_structure, self.ligand_code + '.pdb')

    def help_run(self):
        lpc = LPC(self.native_lig_sdf, self.native_prt_pdb)
        lpc_result = lpc.runLPC()

        pdb = os.path.join(self.dirname,
                           os.path.split(self.native_prt_pdb)[-1])
        shutil.copy(self.native_prt_pdb, pdb)
        # TODO: IOError: openf______E> Cannot open file 3f9mA.seq: No such file or directory
        genProteinEnsemble(self.dirname, pdb, self.ligand_code, lpc_result)

    def run(self):
        self.help_run()
        pass


def test():
    task = PrepareProteinEns("3f9mA")
    luigi.build([task], local_scheduler=True)


def main():
    pass


if __name__ == '__main__':
    import sys
    luigi.build([PrepareProteinEns(sys.argv[1])], local_scheduler=True)
