from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

import os
import shutil
import subprocess


def readContacts(lpc_result):
    result_lines = file(lpc_result).readlines()
    pattern = "Residue      Dist    Surf    HB    Arom    Phob    DC"
    pattern_line_num = -1
    for idx, line in enumerate(result_lines):
        if pattern in line:
            pattern_line_num = idx
    if pattern_line_num == -1: raise "Cannot find contacts in the LPC result"

    # import ipdb; ipdb.set_trace()
    contacts = []
    for idx, line in enumerate(result_lines):
        if idx > pattern_line_num + 1:
            if '----------' in line: break
            contacts.append(line)

    return contacts


# mark the binding residues in the structure sequence
def maskBindingRes(contacts, fasta_seq):
    binding_res_nums = [int(contact.split()[0][0:-1]) for contact in contacts]
    fasta_seq = list(''.join(fasta_seq.split("\n")))
    for res_num in binding_res_nums:
        fasta_seq[res_num - 1] = '-'
    fasta_seq = ''.join(fasta_seq) + "\n"

    return fasta_seq

def isContactRes(contacts, pdb_line):
    binding_res_nums = [int(contact.split()[0][0:-1]) for contact in contacts]
    my_res_num = int(pdb_line.split()[5])
    if my_res_num in binding_res_nums:
        return True
    else:
        return False
    

def genProteinEnsemble(work_dir, prt_code, EdudPrtRmsd = '/home/jaydy/work/working/EdudPrtRmsd'):
    prt_pdb = os.path.abspath("%s/%s.pdb" % (work_dir, prt_code))

    prt_pdb_bk = prt_pdb + '_bk'
    shutil.copyfile(prt_pdb, prt_pdb_bk)  # back up the original protein pdb

    ali_fn = prt_pdb.split('.')[0] + '.ali'

    temp_code = prt_code + '.model'

    lpc_result = EdudPrtRmsd + '/' + prt_code + '/RES1'
    lpc_result = os.path.normpath(lpc_result)

    ################################################################################
    # read the contacts lines from the lpc result file
    contacts = readContacts(lpc_result)

    ################################################################################
    # run babel to convert to fasta
    babel = os.environ['GEAUXDOCK_BABEL']
    if not os.path.exists(babel):
        raise "GEAUXDOCK_BABEL is not properly set"

    cmd = [babel, '-ipdb', prt_pdb, '-ofasta', '-']
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    fasta_seq = "\n".join(out.split("\n")[1:-1]) + "*\n"

    # remove all except for the Ca atoms in the protein
    ca_lines = [line for line in file(prt_pdb_bk) if "CA" in line]
    with open(prt_pdb, 'w') as f:
        my_res_num = 0
        for idx, line in enumerate(ca_lines):
            if not isContactRes(contacts, line):
                my_res_num += 1
                renumbered = "%4s" % str(my_res_num)
                left_part = line[:22]
                right_part = line[26:]
                my_line = left_part + renumbered + right_part
                # print renumbered
                # pass
                f.write(my_line) 


    ################################################################################
    # write to alignment file
    # for PIR format in Modeller, see https://salilab.org/modeller/9v8/manual/node454.html

    ali_lines = []

    ################################################################################
    # header for template
    first_res = ca_lines[0].split()[5]
    last_line = ca_lines[-1]
    chain_id = last_line.split()[4]
    last_res = last_line.split()[5]

    header = ">P1;%s\n%s:%s:%s:%s:%s:%s::::\n" % (temp_code, 'sequence',
                                                  temp_code,
                                                  first_res, chain_id,
                                                  last_res, chain_id)
    ali_lines.append(header + fasta_seq)

    ali_lines.append("\n")

    ################################################################################
    # header for structure
    str_ca_lines = [line for line in file(prt_pdb) if 'CA' in line]
    first_res = str_ca_lines[0].split()[5]
    last_line = str_ca_lines[-1]
    chain_id = last_line.split()[4]
    last_res = last_line.split()[5]
        
    fasta_seq = maskBindingRes(contacts, fasta_seq)  # mask the structure sequence

    header = ">P1;%s\n%s:%s:%s:%s:%s:%s::::\n" % (prt_code, 'structureX',
                                                  prt_code,
                                                  first_res, chain_id,
                                                  last_res, chain_id)
    ali_lines.append(header + fasta_seq)

    ali_fn = os.path.normpath(ali_fn)
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
                  sequence = temp_code)              # code of the target
    a.starting_model= 1                 # index of the first model
    a.ending_model  = 10                 # index of the last model
                                        # (determines how many models to calculate)

    # write the modeller ali seq to file
    # If you also want to see HETATM residues, uncomment this line:
    #env.io.hetatm = True
    code = prt_code
    mdl = model(env, file=code)
    aln = alignment(env)
    aln.append_model(mdl, align_codes=code)
    aln.write(file=code+'.seq')

    a.library_schedule = autosched.slow
    a.max_var_iterations = 300
    a.md_level = refine.slow
    a.repeat_optimization = 3
    a.max_molpdf = 1e6

    # a.very_fast()

    a.make()                            # do the actual comparative modeling

    prt_pdb_for_modeller = prt_pdb + '.mod'
    shutil.copyfile(prt_pdb, prt_pdb_for_modeller)
    shutil.copyfile(prt_pdb_bk, prt_pdb)
