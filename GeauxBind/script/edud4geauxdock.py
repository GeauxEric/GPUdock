#!/usr/bin/env python


import os
import shlex
import subprocess32

import luigi


class Path(luigi.Task):
    ligand_code = luigi.Parameter()

    efindsite_edud_crystal = luigi.Parameter(
        default="/work/jaydy/dat/edud/EDUD-efindiste/efindsite-EDUD-crystal")
    decoy_ligands = luigi.Parameter(
        default="/work/jaydy/dat/edud/decoy-ligands/decoy-ligands/")
    crystal_structure = luigi.Parameter(
        default="/work/jaydy/dat/edud/crystal-structure/")
    crystal_ligands_sdf = luigi.Parameter(
        default="/work/jaydy/dat/edud/crystal-ligands-sdf/")
    active_ligands = luigi.Parameter(
        default="/work/jaydy/dat/edud/active-ligands/")
    work_headquater = luigi.Parameter(
        default="/work/jaydy/working/edud_geaux_inputs/")

    @property
    def work_dir(self):
        path = os.path.join(self.work_headquater, self.ligand_code)
        try:
            os.makedirs(path)
        except:
            pass
        return path


class Extract(Path):
    subset = luigi.Parameter(default='decoy')

    def output(self):
        path = os.path.join(self.subset_work_dir, 'untarred')
        return luigi.LocalTarget(path)

    @property
    def prt_code(self):
        return self.ligand_code[:4]

    @property
    def subset_work_dir(self):
        path = os.path.join(self.work_dir, self.subset)
        try:
            os.makedirs(path)
        except:
            pass
        return path

    def __untar(self, overwrite=True):
        try:
            os.makedirs(self.output().path)
        except:
            pass
        # decoy-ligands
        targz = os.path.join(self.decoy_ligands,
                             self.prt_code + '.tar.gz')
        args = shlex.split(
            '''tar -xf %s -C %s''' % (
                targz, self.output().path
            )
        )
        subprocess32.call(args)

        # efindsite-EDUD-crystal
        targz = os.path.join(self.efindsite_edud_crystal,
                             self.ligand_code,
                             self.ligand_code + '-efindsite.tar.gz')
        args = shlex.split(
            '''tar -xf %s -C %s''' % (
                targz, self.output().path
            )
        )
        subprocess32.call(args)

    def run(self):
        self.__untar()


def test():
    luigi.build(
        [
            Extract("1b9vA")
        ], local_scheduler=True,
    )


def main():
    test()
    pass

if __name__ == '__main__':
    main()
