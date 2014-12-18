import prt_ens
import argparse

################################################################################
#
# generate protein ensembles using Modeller
#
# usage: prt_ens.py [-h] [-d DIRECTORY] [-p PROTEIN]
#
################################################################################

parser = argparse.ArgumentParser(description="generate protein ensemble using Modeller")
parser.add_argument("-d", "--directory", type=str,
                    help="directory for data, must contain the native protein pdb file")
parser.add_argument("-p", "--protein", type=str,
                    help="protein code")
args = parser.parse_args()

work_dir = args.directory
prt_code = args.protein

prt_ens.genProteinEnsemble(work_dir, prt_code)
