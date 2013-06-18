from pymol import cmd
def Simple():
    print "hello PyMol"

cmd.extend("Simple", Simple)
