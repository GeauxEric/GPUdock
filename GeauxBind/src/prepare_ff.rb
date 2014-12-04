#!/home/jaydy/.local/bin/ruby

require 'pry'
require 'optparse'
require_relative './ff'


raise "GEAUXDOCK_PKCOMBU is not properly set\n" unless File.exist? ENV['GEAUXDOCK_PKCOMBU']
raise "GEAUXDOCK_BABEL is not properly set\n" unless File.exist? ENV['GEAUXDOCK_BABEL']
raise "GEAUXDOCK_FF is not properly set\n" unless File.exist? ENV['GEAUXDOCK_FF']

pkcombu = ENV['GEAUXDOCK_PKCOMBU']
babel = ENV['GEAUXDOCK_BABEL']
paramsff = ENV['GEAUXDOCK_FF']

options = {}
OptionParser.new do |opts| 
  opts.banner = "Usage: ruby prepare_ff.rb [options]"

  opts.on("-l LIG", "input ligands in sdf format") do |l| 
    options[:sdf] = l
  end

  opts.on("-i MOLID", "key word for the molecule id, for example: MOLID") do |foo| 
    options[:molid] = foo
  end

  opts.on("-o OUTPUT", "output force field file") do |foo| 
    options[:ff] = foo
  end

  opts.on("-s LIGANDS", "native ligands in sdf format") do |foo| 
    options[:lig] = foo
  end

  opts.on("-a ALIGNMENTS", "alignments data") do |foo| 
    options[:ali] = foo
  end

  opts.on("-p POCKETS", "pockets data") do |foo| 
    options[:pkt] = foo
  end

  opts.on("-t TEMPLATES", "templates proteins in pdb format") do |foo| 
    options[:tmplt] = foo
  end

  opts.on("-n POCKET_NUM", "pocket number to use", Integer) do |foo| 
    options[:pkt_num] = foo
  end

end.parse!



fsdf1 = options[:sdf]
fkey1 = options[:molid]
fout1 = options[:ff]
flig1 = options[:lig]
fali1 = options[:ali]
fpkt1 = options[:pkt]
ftpl1 = options[:tmplt]
fnum1 = options[:pkt_num]


puts "------------------------------------------------------------\n";
puts "                         GeauxDock\n";
puts "                        version 1.0\n\n";
puts "   GPU-accelerated mixed-resolution ligand docking using\n";
puts "                Replica Exchange Monte Carlo\n";
puts "------------------------------------------------------------\n\n";



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
