#!/usr/bin/env python

import shlex
import os
import luigi
import subprocess32
from edud4geauxdock import PrepareLigStateEnsemble


class PrepareLigStageChemicalProps(PrepareLigStateEnsemble):
    esimdock_sdf = luigi.Parameter(
        default=
        "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_sdf")

    def requires(self):
        return PrepareLigStateEnsemble(ligand_code=self.ligand_code,
                                       subset=self.subset)

    def output(self):
        dirname = os.path.dirname(self.requires().output().path)
        path = os.path.join(dirname, self.ligand_code + '_3.sdf')
        return luigi.LocalTarget(path)

    def run(self):
        ifn = self.requires().output().path
        ofn = self.output().path
        cmds = shlex.split('''perl %s -s %s -o %s -i MOLID -c''' %
                           (self.esimdock_sdf, ifn, ofn))
        subprocess32.call(cmds)


class PrepareLigStageForceField(PrepareLigStageChemicalProps):
    prepare_ff = luigi.Parameter(
        default=
        "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/prepare_ff.rb")

    def requires(self):
        return PrepareLigStageChemicalProps(ligand_code=self.ligand_code,
                                            subset=self.subset)

    def output(self):
        path = self.requires().output().path
        out_path = os.path.splitext(path)[0] + '.ff'
        return luigi.LocalTarget(out_path)

    def run(self):
        ifn = self.requires().output().path
        ofn = self.output().path
        untarred_dir = os.path.join(os.path.dirname(ifn), 'untarred')
        ligands_sdf = os.path.join(untarred_dir,
                                   "%s.ligands.sdf" % self.ligand_code)
        align_dat = os.path.join(untarred_dir,
                                 "%s.alignments.dat" % self.ligand_code)
        pockets_dat = os.path.join(untarred_dir,
                                   "%s.pockets.dat" % self.ligand_code)
        templates_pdb = os.path.join(untarred_dir,
                                     "%s.templates.pdb" % self.ligand_code)

        cmds = shlex.split("ruby %s -l %s -i MOLID \
        -o %s -s %s -a %s -p %s -t %s -n 1" % (self.prepare_ff,
                                               ifn,
                                               ofn,
                                               ligands_sdf,
                                               align_dat,
                                               pockets_dat,
                                               templates_pdb, ))

        stdout = subprocess32.check_output(cmds)
        with open(ofn + ".txt", 'w') as ofs:
            ofs.write(stdout)


def main():
    luigi.build([PrepareLigStageForceField("1b9vA")], local_scheduler=True)
    pass


if __name__ == '__main__':
    import sys
    luigi.build([PrepareLigStageForceField(sys.argv[1])], local_scheduler=True)
