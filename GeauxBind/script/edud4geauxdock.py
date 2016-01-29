#!/usr/bin/env python


import os
import logging
import pybel
import shlex
import subprocess32
import luigi

from glob import glob

logging.basicConfig(level=logging.INFO)


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


class PrepareLig1stStage(Extract):
    esimdock_ens = luigi.Parameter(
        default="/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_ens")

    def output(self):
        self.ens_ofn = os.path.join(
            self.subset_work_dir, self.ligand_code + '_2.sdf')
        return luigi.LocalTarget(self.ens_ofn)

    def requires(self):
        return Extract(
            subset=self.subset,
            ligand_code=self.ligand_code
        )

    def aggregate(self):
        """remove hydrogen and add <MOLID>
        """
        untarred_dir = self.requires().output().path
        zincs = glob(os.path.join(untarred_dir, "ZINC*"))
        mols = []
        for zinc in zincs:
            try:
                mols.append(list(pybel.readfile('sdf', zinc)))
            except Exception:
                logging.info("Fail to load zinc ligand %s" % zinc)
        mols = [mol for sub in mols for mol in sub]
        ofn = os.path.join(
            self.subset_work_dir, self.ligand_code + '_1.sdf')
        self.aggregated_ofn = ofn
        ofs = pybel.Outputfile('sdf', ofn, overwrite=True)
        try:
            for mol in mols:
                mol.removeh()
                mol.data['MOLID'] = mol.title
                ofs.write(mol)
        except Exception as detail:
            logging.info(detail)
        finally:
            ofs.close()

    def ens(self):
        cmds = shlex.split(
            '''perl %s -s %s -o %s -i MOLID -c''' % (
                self.esimdock_ens, self.aggregated_ofn, self.ens_ofn))
        subprocess32.call(cmds)

    def run(self):
        self.aggregate()
        self.ens()


def test():
    luigi.build(
        [
            Extract("1b9vA"),
            PrepareLig1stStage("1b9vA"),
        ], local_scheduler=True,
    )


def main():
    test()
    pass

if __name__ == '__main__':
    main()
