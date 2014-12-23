require 'pry'
require 'tempfile'
require 'open3'

# concat the ligand pdb to the end of prt native pdb
def merge(prt_pdb, lig_pdb, complex_pdb)
  prt_lines = File.readlines(prt_pdb)
  lig_lines = File.readlines(lig_pdb)

  prt_atm_lines = prt_lines.select { |line| /ATOM/ =~ line }
  lig_atm_lines = lig_lines.select { |line| /HETATM/ =~ line}

  f = File.open(complex_pdb, 'w')

  f.write("MODEL 1\n")
  prt_atm_lines.each do |line| f.write(line) end
  f.write("TER\n")
  lig_atm_lines.each do |line| f.write(line) end
  f.write("END\n")

  f.close()

end

def runLPC(work_dir, complex_pdb, lpc_bin="/home/jaydy/local/LPC/lpcEx")
  result = File.join(work_dir, 'RES1') 
  File.unlink(result) if File.exist?(result)
  Dir.chdir(work_dir)
  cmd = "#{lpc_bin} 1 #{complex_pdb}"
  stdout_str, stderr_str, status = Open3.capture3(cmd)
  STDERR.puts "Error running #{cmd}\n" unless status.success?
end


# get the protein-ligand contact lines from the LPC result
def readContacts(lpc_result)
  result_lines = File.readlines(lpc_result)
  pattern = "Residue      Dist    Surf    HB    Arom    Phob    DC"
  pattern_line_num = -1
  result_lines.each_with_index do |line, idx| 
    if line.include? pattern
      pattern_line_num = idx
    end
  end
  raise "Cannot find contacts in the LPC result" if pattern_line_num == -1

  contacts = []
  result_lines.each_with_index do |line, idx| 
    if idx > pattern_line_num + 1
      break if line =~ /----------/
      contacts.push(line)
    end
  end

  contacts
end

def zones4Profit(contacts)
  zones = []
  contacts.each do |contact| 
    res = contact.split[0][0...-1]
    zones.push("ZONE #{res}-#{res}\n")
  end
  zones.join
end

def lpcContact(prt_pdb, prt_contact_pdb, lpc_result)

  def contactLigand?(residue, contacts)
    contacts.any { |contact| contact.include? residue }
  end

  def convert(line)
    elems = line.split

    aa_name = elems[3]
    chain_id = elems[4]
    aa_num = elems[5]

    "#{aa_num}#{chain_id}  #{aa_name}"
  end

  contacts = readContacts(lpc_result)

  prt_lines = File.readlines(prt_pdb)
  prt_atm_lines = prt_lines.select { |line| /ATOM/ =~ line }

  contact_lines = prt_atm_lines.select { |line| contacts.any? { |contact| contact.include? convert(line) } }

  f = File.open(prt_contact_pdb, 'w')
  contact_lines.each do |line| f.write(line) end
  f.write("TER\nEND\n")
  f.close()
end


def runProfit(ref_pdb, mob_pdb, zones, profit="~/local/ProFitV3.1/src/profit")
  cmd = "#{profit} -f #{zones} #{ref_pdb} #{mob_pdb}"
  stdout_str, stderr_str, status = Open3.capture3(cmd)
  # binding.pry
  STDERR.puts "Error running #{cmd}\n" unless status.success?

  rms = stdout_str.split("\n").select { |line| line.include? "RMS" }[0].split()[-1]

  rms
end
