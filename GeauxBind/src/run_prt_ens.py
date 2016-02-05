import prt_ens
import argparse
import os

##########################################################################
# usage: run_prt_ens.py [-h] [-d DIRECTORY] [-p PROTEIN] [-c PROTEIN_CODE]
#                       [-l LPC]

# generate protein ensemble using Modeller

# optional arguments:
#   -h, --help            show this help message and exit
#   -d DIRECTORY, --directory DIRECTORY
#                         directory for data, must contain the native protein
#                         pdb file
#   -p PROTEIN, --protein PROTEIN
#                         protein pdb file
#   -c PROTEIN_CODE, --protein_code PROTEIN_CODE
#                         protein code
#   -l LPC, --lpc LPC     lpc result
##########################################################################

parser = argparse.ArgumentParser(
    description="generate protein ensemble using Modeller")
parser.add_argument(
    "-d",
    "--directory",
    type=str,
    help="directory for data, must contain the native protein pdb file")
parser.add_argument("-p", "--protein", type=str, help="protein pdb file")
parser.add_argument("-c", "--protein_code", type=str, help="protein code")
parser.add_argument("-l", "--lpc", type=str, help="lpc result")

args = parser.parse_args()

work_dir = args.directory
prt_pdb = args.protein
prt_code = args.protein_code
lpc_result = args.lpc

os.chdir(work_dir)
prt_ens.genProteinEnsemble(work_dir, prt_pdb, prt_code, lpc_result)
