from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

import os
import shutil
import subprocess

def genProteinEnsemble(work_dir, prt_code):
    prt_pdb = "%s/%s.pdb" % (work_dir, prt_code)
    prt_pdb_bk = prt_pdb + '_bk'
    ali_fn = prt_pdb.split('.')[0] + '.ali'

    ################################################################################
    # run babel to convert to fasta
    babel = os.environ['GEAUXDOCK_BABEL']
    if not os.path.exists(babel):
        raise "GEAUXDOCK_BABEL is not properly set"

    cmd = [babel, '-ipdb', prt_pdb, '-ofasta', '-']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    fasta_seq = "\n".join(out.split("\n")[1:-1]) + "*"

    # remove all except for the Ca atoms in the protein
    shutil.copyfile(prt_pdb, prt_pdb_bk)  # back up the original protein pdb
    ca_lines = [line for line in file(prt_pdb_bk) if "CA" in line]
    with open(prt_pdb, 'w') as f:
        for line in ca_lines:
            f.write(line)

    ################################################################################
    # write to alignment file
    ali_lines = []
    ali_lines.append(">P1;model\nsequence:model:::::::0.00:0.00\n" + fasta_seq)
    ali_lines.append("\n\n")
    header = ">P1;%s\nstructureX:%s.pdb:1:A:::undefined:undefined:-1.00:-1.00\n" % (prt_code, prt_code)
    ali_lines.append(header + fasta_seq)

    with open(ali_fn, 'w') as f:
        for line in ali_lines:
            f.write(line)

    ################################################################################
    # modeller
    env = environ()  # create a new MODELLER environment to build this model in

    # directories for input atom files
    env.io.atom_files_directory = [work_dir]

    a = automodel(env,
                  alnfile  = ali_fn,     # alignment filename
                  knowns   = prt_code,              # codes of the templates
                  sequence = prt_code)              # code of the target
    a.starting_model= 1                 # index of the first model
    a.ending_model  = 10                 # index of the last model
                                        # (determines how many models to calculate)

    a.library_schedule = autosched.slow
    a.max_var_iterations = 300
    a.md_level = refine.slow
    a.repeat_optimization = 3
    a.max_molpdf = 1e6

    # a.very_fast()

    a.make()                            # do the actual comparative modeling

    shutil.copyfile(prt_pdb_bk, prt_pdb)
