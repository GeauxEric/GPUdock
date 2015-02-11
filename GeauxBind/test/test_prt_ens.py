import unittest

import subprocess
import prt_ens
import os

class TestProteinEnsemble(unittest.TestCase):

    def setUp(self):
        print "run the test_lpc.rb test for the LPC result of docking 1b9vA"
        print "--------------------------------------------------------------------------------"
        cmd = ['ruby', 'test_lpc.rb']
        subprocess.call(cmd)

    def test_shuffle(self):
        # generate ensembles for 1b9vA
        work_dir = "../data"
        prt_pdb = "../data/1b9vA.pdb"
        prt_code = "1b9vA"
        lpc_result = "../data/RES1"
        os.chdir(work_dir)
        prt_ens.genProteinEnsemble(work_dir, prt_pdb, prt_code, lpc_result)

if __name__ == '__main__':
    unittest.main()
