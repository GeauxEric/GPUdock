require 'test/unit'
require_relative 'prepare_sdf'

class PrepareSdfTest < Test::Unit::TestCase

  def test_removeH
    ifn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158"
    ofn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_1.sdf"
    babel = "/home/jaydy/local/bin/babel"

    removeH(ifn, ofn, babel=babel)
  end

  def test_addMolid
    ifn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_1.sdf"
    ofn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_2.sdf"
    addMolid(ifn, ofn)
  end

  def test_run_sdf
    ifn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_2.sdf"
    ofn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_3.sdf"
    perl_script = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_sdf"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -c"
    sdf_runs_well = system cmd
    raise "failure in running #{perl_script}" unless sdf_runs_well
  end

  def test_run_ens
    ifn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_3.sdf"
    ofn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158_4.sdf"
    perl_script = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/src/esimdock_ens"
    cmd = "perl #{perl_script} -s #{ifn} -o #{ofn} -i MOLID -n 50"
    ens_runs_well = system cmd
    raise "failure in running #{perl_script}" unless ens_runs_well
  end

end
