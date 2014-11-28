require 'test/unit'
require 'tempfile'
require_relative 'prepare_sdf'


class PrepareSdfTest < Test::Unit::TestCase

  ################################################################################
  # removeH

  # on single compound
  def test_aa_removeH
    puts "\n################################################################################\nremove H\n"
    babel = "/home/jaydy/local/bin/babel"

    ifn = "../data/ZINC00002158"
    ofn = "../data/ZINC00002158_1.sdf"

    removeH(ifn, ofn, babel=babel)
    puts "write to #{ofn}\n"
  end


  ################################################################################
  # addMolid2File

  # on single compound
  def test_ba_addMolid2File
    puts "\n################################################################################\nadd MOLID\n"
    ifn = "../data/ZINC00002158_1.sdf"
    ofn = "../data/ZINC00002158_2.sdf"
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf

  # on single compound
  def test_ca_run_sdf
    puts "\n################################################################################\nrun esimdock_sdf\n"
    ifn = "../data/ZINC00002158_2.sdf"
    ofn = "../data/ZINC00002158_3.sdf"
    perl_script = "./esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    puts "\n" + cmd + "\n"

    runs_well = system cmd
    if runs_well
      puts "write to #{ofn}"
    else
      raise "failure in running #{perl_script}\n"
    end
  end


  ################################################################################
  # run esimdock_ens

  # on single compound
  def test_da_run_ens
    puts "\n################################################################################\nrun esimdock_ens\n"
    ifn = "../data/ZINC00002158_3.sdf"
    ofn = "../data/ZINC00002158_4.sdf"
    perl_script = "./esimdock_ens"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -n 50"
    puts "\n" + cmd + "\n"

    runs_well = system cmd
    if runs_well
      puts "write to #{ofn}"
    else
      raise "failure in running #{perl_script}\n"
    end
  end

  ################################################################################
  # run prepare_ff

  # on single compound
  def test_ea_run_prepare_ff
    puts "\n################################################################################\nrun prepare_ff\n"
    ifn = "../data/ZINC00002158_4.sdf"
    ofn = "../data/ZINC00002158_4.ff"

    perl_script = "./prepare_ff"
    cmd = "perl #{perl_script} -l #{ifn} -i MOLID -o #{ofn} \
-s ../data/1b9vA.ligands.sdf -a ../data/1b9vA.alignments.dat \
-p ../data/1b9vA.pockets.dat -t ../data/1b9vA.templates.pdb -n 1"

    puts "\n" + cmd + "\n"

    runs_well = system cmd
    if runs_well
      puts "write to #{ofn}"
    else
      raise "failure in running \n #{cmd}\n"
    end
  end

end
