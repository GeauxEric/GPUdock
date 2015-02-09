require 'test/unit'
require 'tempfile'
require 'open3'
require_relative '../src/prepare_sdf'


class PrepareSdfTest < Test::Unit::TestCase

  ################################################################################
  # removeH

  # on multiple compound
  def test_ab_removeH
    babel = "/home/jaydy/local/bin/babel"

    ifn1 = "../data/ZINC00002158"
    ifn2 = "../data/ZINC00089285"
    ifn3 = "../data/ZINC01572843"
    ifns = [ifn1, ifn2, ifn3]

    tmp = Tempfile.new("/ZINC")
    tmp_path = tmp.path
    tmp.close

    concatFiles(ifns, tmp_path)

    ofn = "../data/ZINC_1.sdf"
    removeH(tmp_path, ofn, babel=babel)

    tmp.unlink
  end

  ################################################################################
  # addMolid2File

  # on multiple compounds
  def test_bb_addMolid2File
    ifn = "../data/ZINC_1.sdf"
    ofn = "../data/ZINC_2.sdf"
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf

  # on multiple compounds
  def test_cb_run_sdf
    ifn = "../data/ZINC_2.sdf"
    ofn = "../data/ZINC_3.sdf"
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
    ifn = "../data/ZINC_3.sdf"
    ofn = "../data/ZINC_4.sdf"
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

  
  ################################################################################
  # run the ruby version of prepare_ff
  def test_fb_run_prepare_ff
    ifn = "../data/ZINC_4.sdf"
    ofn = "../data/ZINC_5.ff"

    ruby_script = "../src/prepare_ff.rb"
    cmd = "ruby #{ruby_script} -l #{ifn} -i MOLID -o #{ofn} \
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

end
