#!/home/jaydy/.local/bin/ruby

require 'pry'
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
kde = prepareKDE(eligible_ligs, babel=babel)

puts "done\n\n";
 
puts "Calculating pocket-specific potential ... \n";
psp = preparePSP(fali1, paramsff, fpkt1, ftpl1, eligible_ligs, fnum1, kde, babel=babel)
puts "done\n\n"

print "Calculating position restraints for docking compounds:\n\n";
mcs = prepareMCS(fsdf1, eligible_ligs, fkey1, pkcombu=pkcombu)
print "done\n"

output_lines = kde.join("\n") + "\n" + psp.join("\n") + "\n" + mcs.join("\n")

File.open(fout1, 'w') do |file| 
  file.write(output_lines)
end

ref_ifn = "../data/ZINC_single_correct.ff"
ref_lines = File.open(ref_ifn).readlines.join

raise "Wroing generating .ff" unless output_lines == ref_lines
