require 'test/unit'
require 'tempfile'
require 'open3'
require_relative '../src/prepare_sdf'


class PrepareSdfTest < Test::Unit::TestCase

  ################################################################################
  # removeH

  # on single compound
  def test_aa_removeH
    babel = "/home/jaydy/local/bin/babel"

    ifn = "../data/ZINC00002158"
    ofn = "../data/ZINC00002158_1.sdf"

    cmd = "#{babel} -isdf #{ifn} -osdf #{ofn} -d"
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
  # addMolid2File

  # on single compound
  def test_ba_addMolid2File
    ifn = "../data/ZINC00002158_1.sdf"
    ofn = "../data/ZINC00002158_2.sdf"
    addMolid2File(ifn, ofn)
    puts "\nwrite to\t#{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf

  # on single compound
  def test_ca_run_sdf
    ifn = "../data/ZINC00002158_2.sdf"
    ofn = "../data/ZINC00002158_3.sdf"
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

  # on single compound
  def test_da_run_ens
    ifn = "../data/ZINC00002158_3.sdf"
    ofn = "../data/ZINC00002158_4.sdf"
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
  # run prepare_ff

  # on single compound
  def test_ea_run_prepare_ff
    ifn = "../data/ZINC00002158_4.sdf"
    ofn = "../data/ZINC00002158_4.ff"
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

end
