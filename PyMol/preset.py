from pymol import cmd

# set the background to white
def white_bg():
    set bg_rgb=[1,1,1]      # these are [red, blue, green] components

cmd.extend("white_bg", white_bg)

