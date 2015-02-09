require 'open3'
require_relative 'prepare_sdf'


class EdudForceField
  attr_accessor :work_dir, :prt_id

  ################################################################################
  # removeH
  # on multiple compound
  def run_removeH()
    babel = "/home/jaydy/local/bin/babel"

    chembl_fns = Dir.glob(@work_dir + '/CHEMBL*')

    ifns = chembl_fns.sample(6)
    # ifns = chembl_fns

    tmp = Tempfile.new('/CHEMBL')
    tmp_path = tmp.path
    tmp.close

    concatFiles(ifns, tmp_path)

    ofn = @work_dir + '_1.sdf'
    removeH(tmp_path, ofn, babel=babel)

    tmp.unlink

  end

  ################################################################################
  # addMolid2File
  # on multiple compounds
  def run_addMolid2File()
    ifn = @work_dir + '_1.sdf'
    ofn = @work_dir + '_2.sdf'
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf
  # on multiple compounds
  def run_sdf(perl_script)
    ifn = @work_dir + '_2.sdf'
    ofn = @work_dir + '_3.sdf'
    # perl_script = "../src/esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    puts "\nRunning\t\t#{cmd}\n"
    stdout_str, stderr_str, status = Open3.capture3(cmd)

    if status.success?
      puts "write to\t#{ofn}\n"
    else
      STDERR.puts "Error running #{cmd}\n"
      exit 1
    end
  end

  ################################################################################
  # run esimdock_ens
  def run_ens(perl_script)
    ifn = @work_dir + '_3.sdf'
    ofn = @work_dir + '_4.sdf'
    # perl_script = "../src/esimdock_ens"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -n 50"
    puts "\nRunning\t\t#{cmd}\n"
    stdout_str, stderr_str, status = Open3.capture3(cmd)

    if status.success?
      puts "write to\t#{ofn}\n"
    else
      STDERR.puts "Error running #{cmd}\n"
      exit 1
    end
  end


  
  ################################################################################
  # run the ruby version of prepare_ff
  def run_prepare_ff(ruby_script)
    ifn = @work_dir + '_4.sdf'
    ofn = @work_dir + '_4.ff'

    ligands_sdf = File.join(@work_dir, @prt_id + '.ligands.sdf')
    ali = File.join(@work_dir, @prt_id + '.alignments.dat')
    pocket = File.join(@work_dir, @prt_id + '.pockets.dat')
    templates = File.join(@work_dir, @prt_id + '.templates.pdb')

    cmd = <<ff
ruby #{ruby_script} -l #{ifn} -i MOLID -o #{ofn} -s #{ligands_sdf} -a #{ali} -p #{pocket} -t #{templates} -n 1
ff

    puts "\nRunning\t\t#{cmd}\n"
    stdout_str, stderr_str, status = Open3.capture3(cmd)
    if status.success?
      puts "write to\t#{ofn}\n"
    else
      STDERR.puts "Error running #{cmd}\n"
      exit 1
    end
  end

end
