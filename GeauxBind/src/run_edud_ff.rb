require 'optparse'
require_relative 'edud_ff'

options = {}
OptionParser.new do |opts| 
  opts.banner = "Usage: ruby run_edud_ff.rb [options]"

  opts.on("-l LIG", "input ligand id") do |l| 
    options[:lig] = l
  end

  opts.on("-p PRT", "input protein id") do |p| 
    options[:prt] = p
  end

  opts.on("-n POCKET", Integer, "pocket number") do |n| 
    options[:pocket] = n
  end

  # list of the paths of the anchor structures
  opts.on("-a ", "anchor structures") do |a| 
    options[:anchors] = a
  end

end.parse!

# lig_id = '1bcd'
# prt_id = '1bcdA'

lig_id = options[:lig]
prt_id = options[:prt]
pocket_num = options[:pocket]
anchors_ifns = options[:anchors]


puts "------------------------------------------------------------\n";
puts "                    GeauxDock Force-Field\n";
puts "                        version 0.1\n\n";
puts "------------------------------------------------------------\n\n";

################################################################################
# prepare the files
edud_dat = '/work/jaydy/dat/EDUD'
efindsite_edud = '/work/jaydy/dat/EDUD/efindsite-EDUD-crystal/'

env = '/home/jaydy/work/working/EdudCrystalLig/'
esimdock_sdf = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_sdf'
esimdock_ens = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_ens'
ruby_prepare_ff = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/prepare_ff.rb'

lig_sdf = File.join(edud_dat, lig_id + '.sdf')
prt_pdb = File.join(edud_dat, prt_id + '.pdb')
edud_tar = File.join(edud_dat, lig_id + '.tar.gz')
efindsite_tar = File.join(efindsite_edud, prt_id, prt_id + '-efindsite.tar.gz')

work_dir = File.join(env, lig_id)
FileUtils.mkdir_p(work_dir)
FileUtils.cp(lig_sdf, work_dir)
FileUtils.cp(prt_pdb, work_dir)
`tar xf #{edud_tar} -C #{work_dir}`
`tar xf #{efindsite_tar} -C #{work_dir}`

anchors_ifns = Dir.glob(File.join(work_dir, 'CHEMBL*'))
ifns = anchors_ifns[0, 6]

################################################################################
# run
ff_job = EdudForceField.new()
ff_job.work_dir = work_dir
ff_job.prt_id = prt_id

ff_job.run_removeH(ifns)
ff_job.run_addMolid2File()
ff_job.run_sdf(esimdock_sdf)
ff_job.run_ens(esimdock_ens)
ff_job.run_prepare_ff(ruby_prepare_ff, pocket_num)
