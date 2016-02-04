#!/usr/bin/env python

import luigi

from edud4geauxdock2 import PrepareLigStageForceField
from edud4geauxdock import PrepareLigStateEnsemble


class CheckEns(luigi.Task):
    def output(self):
        path = "./pbs/finished.txt"
        return luigi.LocalTarget(path)

    def requires(self):
        pass

    def run(self):
        names = [line.rstrip() for line in file("./pbs/edud_input_jobs.txt")]
        completes, incompletes = [], []
        for name in names:
            task = PrepareLigStateEnsemble(name)
            if task.complete():
                completes.append(name)
            else:
                incompletes.append(name)
        with open(self.output().path, 'w') as ofs:
            for name in completes:
                ofs.write("{}\n".format(name))


class CheckCompletion(luigi.Task):
    @property
    def jobs(self):
        names = [line.rstrip() for line in file("./pbs/finished.txt")]
        return names

    def requires(self):
        return CheckEns()

    def output(self):
        pass

    def run(self):
        completes, incompletes = [], []
        for name in self.jobs:
            task = PrepareLigStageForceField(name)
            if task.complete():
                completes.append(task)
            else:
                incompletes.append(task)

        print("{} completes and {} incompletes".format(
            len(completes), len(incompletes)))


def main():
    luigi.build([CheckCompletion(), CheckEns()], local_scheduler=True)


if __name__ == '__main__':
    main()
