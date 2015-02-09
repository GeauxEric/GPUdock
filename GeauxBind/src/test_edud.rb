require 'test/unit'
require 'tempfile'
require 'open3'
require_relative 'prepare_sdf'


class PrepareSdfTest < Test::Unit::TestCase
  @@lig = '1b9vA'
  @@prt = '1b9v'
  @@work_dir = '/home/jaydy/work/working/EdudCrystalLig/'
  @@env = File.join(@@work_dir, @@prt)

  ################################################################################
  # removeH

  # on multiple compound
  def test_ab_removeH
    babel = "/home/jaydy/local/bin/babel"

    chembl_fns = Dir.glob(@@env + '/CHEMBL*')

    ifns = chembl_fns.sample(6)

    tmp = Tempfile.new('/CHEMBL')
    tmp_path = tmp.path
    tmp.close

    concatFiles(ifns, tmp_path)

    ofn = @@env + '_1.sdf'
    removeH(tmp_path, ofn, babel=babel)

    tmp.unlink
  end

  ################################################################################
  # addMolid2File

  # on multiple compounds
  def test_bb_addMolid2File
    ifn = @@env + '_1.sdf'
    ofn = @@env + '_2.sdf'
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf

  # on multiple compounds
  def test_cb_run_sdf
    ifn = @@env + '_2.sdf'
    ofn = @@env + '_3.sdf'
    perl_script = "../src/esimdock_sdf"
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

  # on multiple compounds
  def test_db_run_ens
    ifn = @@env + '_3.sdf'
    ofn = @@env + '_4.sdf'
    perl_script = "../src/esimdock_ens"
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


=begin
  ################################################################################
  # run the perl version of prepare_ff

  def test_eb_run_prepare_ff
    ifn = "../data/ZINC_4.sdf"
    ofn = "../data/ZINC_4.ff"

    perl_script = "../src/prepare_ff"
    cmd = "perl #{perl_script} -l #{ifn} -i MOLID -o #{ofn} \
-s ../data/1b9vA.ligands.sdf -a ../data/1b9vA.alignments.dat \
-p ../data/1b9vA.pockets.dat -t ../data/1b9vA.templates.pdb -n 1"

    puts "\nRunning\t\t#{cmd}\n"
    stdout_str, stderr_str, status = Open3.capture3(cmd)
    if status.success?
      puts "write to\t#{ofn}\n"
    else
      STDERR.puts "Error running #{cmd}\n"
      exit 1
    end
  end
=end
  
  ################################################################################
  # run the ruby version of prepare_ff
  def test_fb_run_prepare_ff
    ifn = @@env + '_4.sdf'
    ofn = @@env + '_4.ff'

    ligands_sdf = File.join(@@env, @@lig + '.ligands.sdf')
    ali = File.join(@@env, @@lig + '.alignments.dat')
    pocket = File.join(@@env, @@lig + '.pockets.dat')
    templates = File.join(@@env, @@lig + '.templates.pdb')

    ruby_script = "../src/prepare_ff.rb"
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
