require_relative 'edud_ff'

edud_dat = '/work/jaydy/dat/EDUD'
efindsite_edud = '/work/jaydy/dat/EDUD/efindsite-EDUD-crystal/'

env = '/home/jaydy/work/working/EdudCrystalLig/'
esimdock_sdf = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_sdf'
esimdock_ens = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_ens'
ruby_prepare_ff = '/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/prepare_ff.rb'

lig_id = '1bcd'
prt_id = '1bcdA'


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

ff_job = EdudForceField.new()
ff_job.work_dir = work_dir
ff_job.prt_id = prt_id

ff_job.run_removeH()
ff_job.run_addMolid2File()
ff_job.run_sdf(esimdock_sdf)
ff_job.run_ens(esimdock_ens)
ff_job.run_prepare_ff(ruby_prepare_ff)

