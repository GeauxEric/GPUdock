require 'test/unit'
require_relative 'prepare_sdf'

class PrepareSdfTest < Test::Unit::TestCase

  def test_removeH
    ifn = "../data/ZINC00002158"
    ofn = "../data/ZINC00002158_1.sdf"
    babel = "/home/jaydy/local/bin/babel"

    removeH(ifn, ofn, babel=babel)
  end

  def test_addMolid2File
    ifn = "../data/ZINC00002158_1.sdf"
    ofn = "../data/ZINC00002158_2.sdf"
    addMolid2File(ifn, ofn)
  end

  def test_run_sdf
    ifn = "../data/ZINC00002158_2.sdf"
    ofn = "../data/ZINC00002158_3.sdf"
    perl_script = "./esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    sdf_runs_well = system cmd
    raise "failure in running #{perl_script}" unless sdf_runs_well
  end

  def test_run_ens
    ifn = "../data/ZINC00002158_3.sdf"
    ofn = "../data/ZINC00002158_4.sdf"
    perl_script = "./esimdock_ens"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -n 50"
    ens_runs_well = system cmd
    raise "failure in running #{perl_script}" unless ens_runs_well
  end

  def test_run_ff
    ifn = "../data/ZINC00002158_4.sdf"
    ofn = "../data/ZINC00002158_4.ff"
    perl_script = "./prepare_ff"
    cmd = "perl #{perl_script} -l #{ifn} \
        -i MOLID \
        -o #{ofn} \
        -s ../data/1b9vA.ligands.sdf \
        -a ../data/1b9vA.alignments.dat \
        -p ../data/1b9vA.pockets.dat \
        -t ../data/1b9vA.templates.pdb \
        -n 1"
    ff_runs_well = system cmd
    raise "failure in running #{perl_script}" unless ff_runs_well
  end

end
