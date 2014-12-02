#!/home/jaydy/.local/bin/ruby

require 'optparse'
require_relative '../src/ff'

puts "------------------------------------------------------------\n";
puts "                         GeauxDock\n";
puts "                        version 1.0\n\n";
puts "   GPU-accelerated mixed-resolution ligand docking using\n";
puts "                Replica Exchange Monte Carlo\n";
puts "------------------------------------------------------------\n\n";


raise "GEAUXDOCK_PKCOMBU is not set\n" unless File.exist? ENV['GEAUXDOCK_PKCOMBU']
raise "GEAUXDOCK_BABEL is not set\n" unless File.exist? ENV['GEAUXDOCK_BABEL']
raise "GEAUXDOCK_FF is not set\n" unless File.exist? ENV['GEAUXDOCK_FF']

pkcombu = ENV['GEAUXDOCK_PKCOMBU']
babel = ENV['GEAUXDOCK_BABEL']
paramsff = ENV['GEAUXDOCK_FF']

# default arguments
# perl ../src/prepare_ff -l ../data/ZINC00002158_4.sdf -i MOLID -o ../data/ZINC00002158_4.ff -s ../data/1b9vA.ligands.sdf -a ../data/1b9vA.alignments.dat -p ../data/1b9vA.pockets.dat -t ../data/1b9vA.templates.pdb -n 1

fsdf1 = '../data/ZINC00002158_4.sdf'
fkey1 = 'MOLID'
fout1 = '../data/ZINC00002158_4.ff'
flig1 = '../data/1b9vA.ligands.sdf'
fali1 = '../data/1b9vA.alignments.dat'
fpkt1 = '../data/1b9vA.pockets.dat'
ftpl1 = '../data/1b9vA.templates.pdb'
fnum1 = 1

puts "Preparing data for KDE ... \n";
eligible_ligs = getEligibleLigs(flig1, fnum1)
kde_pts, kde_atom_num = prepareKDE(eligible_ligs, babel=babel)
raise "error calculating KDE" unless kde_pts[0] == "KDE C.2 33.143 -12.2923 66.0406\n"
puts "done\n\n";
 
puts "Calculating pocket-specific potential ... \n";

psp = preparePSP(fali1, paramsff, fpkt1, ftpl1, eligible_ligs, fnum1, kde_pts, kde_atom_num, babel=babel)
puts kde_pts
puts psp
