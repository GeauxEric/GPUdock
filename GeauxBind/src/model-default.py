# Comparative modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['../data/']

a = automodel(env,
              alnfile  = '../data/align-ligand.ali',     # alignment filename
              knowns   = '1j4hA',              # codes of the templates
              sequence = 'model')              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
                                    # (determines how many models to calculate)

a.library_schedule = autosched.slow
a.max_var_iterations = 300
a.md_level = refine.slow
a.repeat_optimization = 3
a.max_molpdf = 1e6

a.make()                            # do the actual comparative modeling
