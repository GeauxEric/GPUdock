require 'test/unit'
require 'tempfile'
require_relative 'prepare_sdf'


class PrepareSdfTest < Test::Unit::TestCase

  ################################################################################
  # removeH
  def test_a_removeH
    babel = "/home/jaydy/local/bin/babel"

    ifn = "../data/ZINC00002158"
    ofn = "../data/ZINC00002158_1.sdf"

    removeH(ifn, ofn, babel=babel)
  end

  def test_b_removeH
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
    # run babel on a sdf file containing more than one compounds
    removeH(tmp_path, ofn, babel=babel)

    tmp.unlink
  end

  ################################################################################
  # addMolid2File
  def test_c_addMolid2File
    ifn = "../data/ZINC00002158_1.sdf"
    ofn = "../data/ZINC00002158_2.sdf"
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  def test_d_addMolid2File
    ifn = "../data/ZINC_1.sdf"
    ofn = "../data/ZINC_2.sdf"
    addMolid2File(ifn, ofn)
    puts "\nwrite to #{ofn}\n"
  end

  ################################################################################
  # run esimdock_sdf
  def test_e_run_sdf
    ifn = "../data/ZINC00002158_2.sdf"
    ofn = "../data/ZINC00002158_3.sdf"
    perl_script = "./esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    puts "\n" + cmd + "\n"

    sdf_runs_well = system cmd
    if sdf_runs_well
      puts "write to #{ofn}"
    else
      raise "failure in running #{perl_script}\n"
    end
  end

  def test_f_run_sdf
    ifn = "../data/ZINC_2.sdf"
    ofn = "../data/ZINC_3.sdf"
    perl_script = "./esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    puts "\n" + cmd + "\n"

    sdf_runs_well = system cmd
    if sdf_runs_well
      puts "write to #{ofn}"
    else
      raise "failure in running #{perl_script}\n"
    end
  end


end
