require 'test/unit'
require 'tempfile'
require 'open3'
require_relative '../src/lpc'


class PrepareSdfTest < Test::Unit::TestCase

  def test_aa_merge
    prt_pdb = "../data/1b9vA.pdb"
    lig_sdf = "../data/1b9v.sdf"
    lig_pdb = "../data/1b9v.pdb"
    `babel -isdf #{lig_sdf} -opdb #{lig_pdb}`

    complex_pdb = "../data/1b9vA.complex.pdb"
    merge(prt_pdb, lig_pdb, complex_pdb)
  end

  def test_ab_runLPC
    lpcEx = "/home/jaydy/local/LPC/lpcEx"
    complex_pdb = "../data/1b9vA.complex.pdb"
    work_dir = "../data"
    runLPC(work_dir, complex_pdb, lpc_bin=lpcEx)
  end

  def test_ac_readContacts
    work_dir = "../data"
    lpc_result = File.join(work_dir, "RES1")
    readContacts(lpc_result)
  end

  def test_ac_zones4Profit
    work_dir = "../data"
    lpc_result = File.join(work_dir, "RES1")
    contacts = readContacts(lpc_result)
    zones4Profit(contacts)
  end

  def test_ad_lpcContact
    prt_pdb = "../data/1b9vA.pdb"
    prt_contact_pdb = "../data/1b9vA.contact.pdb"
    work_dir = "../data"
    lpc_result = File.join(work_dir, "RES1")
    lpcContact(prt_pdb, prt_contact_pdb, lpc_result)
  end

  def test_ae_runProfit
    work_dir = "../data"
    lpc_result = File.join(work_dir, "RES1")
    contacts = readContacts(lpc_result)
    zones = zones4Profit(contacts)

    script_fn = File.join(work_dir, "RES1.fit")
    f = File.open(script_fn, 'w')
    f.write(zones)
    f.write("FIT\n")
    f.close()
    
  end

  def test_af_runProfit
    ref_pdb = "../data/reference.pdb"
    mob_pdb = "../data/mobile.pdb"
    rms = runProfit(ref_pdb, mob_pdb)
    puts "\nRMSD\t #{rms}"
  end
  

end
