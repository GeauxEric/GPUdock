require 'test/unit'
require_relative 'prepare_sdf'

class PrepareSdfTest < Test::Unit::TestCase

  def test_removeH
    babel = "/home/jaydy/local/bin/babel"

    ifn = "../data/ZINC00002158"
    ofn = "../data/ZINC00002158_1.sdf"

    removeH(ifn, ofn, babel=babel)

    ifn = "../data/ZINC00089285"
    ofn = "../data/ZINC00089285_1.sdf"

    removeH(ifn, ofn, babel=babel)
  end

  def test_addMolid2File
    ifn = "../data/ZINC00002158_1.sdf"
    ofn = "../data/ZINC00002158_2.sdf"
    addMolid2File(ifn, ofn)

    ifn = "../data/ZINC00089285_1.sdf"
    ofn = "../data/ZINC00089285_2.sdf"
    addMolid2File(ifn, ofn)
  end

  def test_addMultipleMolids
    ifns = ["../data/ZINC00002158_1.sdf", "../data/ZINC00089285_1.sdf"]
    ofn = "../data/ZINC0000.sdf"
    File.delete(ofn) if File.exist?(ofn)
    addMultipleMolids(ifns, ofn)
    print "write to #{ofn}\n"
  end

  def test_run_sdf
    ifn = "../data/ZINC0000.sdf"
    ofn = "../data/ZINC0000_1.sdf"
    perl_script = "./esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    puts cmd + "\n"
    sdf_runs_well = system cmd
    raise "failure in running #{perl_script}\n" unless sdf_runs_well
  end


end
