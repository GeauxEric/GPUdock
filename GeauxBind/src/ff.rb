require 'pry'
require 'tempfile'
require 'open3'

# get ligand atom name and coord in the mol2 format
def getLigAtomsFromMol2 (mol2_lines)
  lig_atoms = []
  atom_line_is_found = false
  mol2_lines.each do |line|
    atom_line_is_found = false if line == '@<TRIPOS>BOND'
    if atom_line_is_found
      x = line[16, 10].to_f.to_s
      y = line[26, 10].to_f.to_s
      z = line[36, 10].to_f.to_s
      atom_name = line[47, 7].delete(' ')
      lig_atoms.push([atom_name, x, y, z])
    end
    atom_line_is_found = true if line == '@<TRIPOS>ATOM'
  end
  return lig_atoms
end

# use the seperator_requirement method to seperate the text into several sections
def seperateLines(lines, seperator_requirement)
  splitted = []
  section = []
  lines.each do |line|
    section.push(line)
    if seperator_requirement.call(line)
      splitted.push(section)
      section = []
    end
  end

  splitted
end

def getPocketCenterCoords(fpkt1, fnum1)
  lines = File.readlines(fpkt1).each{ |line| line.chomp! }
  seperator_line_requirement = Proc.new do |line| line == "TER" end 
  pockets = seperateLines(lines, seperator_line_requirement)

  def getPocketNum(pocket)
    num = -1
    pocket.each do |line| 
      num = line.split[1].to_i if line.start_with?("POCKET")
    end
    num
  end

  def getCoords(pocket)
    pocket.each do |line| 
      if line.start_with?("CENTER")
        return line.split[1, 3]
      end
    end
  end

  my_center = []
  pockets.each do |pocket| 
    my_num = getPocketNum(pocket)
    if my_num == fnum1
      my_center = getCoords(pocket)
    end
  end

  if not my_center.empty?
    return 'CENTER ' + my_center.join(' ')
  else
    raise "not center found for" + fpkt1
  end
end

# read a .ligands.sdf file, get the ligands with required pocket number
def getEligibleLigs(flig1, fnum1)

  def getPocketNum(lig)
    pocket = 0
    lig.each_with_index do |line, idx|
      pocket = lig[idx + 1] if line =~ /EFINDSITE_POCKET/
    end

    pocket.to_i
  end

  def getMolid(lig)
    molid = ''
    lig.each_with_index do |line, idx|
      molid = lig[idx + 1] if line =~ /MOLID/
    end

    molid
  end

  def requiredPocketNum?(lig, fnum)
    getPocketNum(lig) == fnum ? true : false
  end

  lig01 = File.readlines(flig1).each{ |line| line.chomp! }

  seperator_line_requirement = Proc.new do |line| line == "$$$$" end

  ligs = seperateLines(lig01, seperator_line_requirement)

  eligible_ligs = Hash.new
  ligs.each do |lig|
    if requiredPocketNum?(lig, fnum1)
      molid = getMolid(lig)
      eligible_ligs[molid] = lig
    end
  end

  return eligible_ligs
end

# convert from sdf to mol2
# grep the ligand coords and name in the mol2 file
def prepareKDE(eligible_ligs, babel='')

  lines_to_sdf = eligible_ligs.values.flatten
  sdf_tmp = Tempfile.new("kde_sdf")
  path = sdf_tmp.path
  sdf_tmp.close
  File.open(path, 'w') do |file|
    file.write(lines_to_sdf.join("\n"))
  end

  # convert from sdf to mol2
  cmd = "#{babel} -isdf #{path} -omol2 -"
  stdout_str, stderr_str, status = Open3.capture3(cmd)
  unless status.success?
    STDERR.puts "Error running #{cmd}\n"
  end

  sdf_tmp.unlink

  mol2 = stdout_str.split("\n")
  lig_atoms = getLigAtomsFromMol2(mol2)

  kde = []
  lig_atoms.each do |atom| 
    n, x, y, z = atom
    kde.push('KDE ' + [n, x, y, z].join(' '))
  end

  return kde
end

def preparePSP(fali1, paramsff, fpkt1, ftpl1, eligible_ligs, fnum1, kde_pts, babel='')
  ################################################################################
  # deal with .alignments.dat
  ali01 = File.open(fali1, 'r').readlines.each do |line| line.chomp! end
  seperator_line_requirement = Proc.new do |line| line == "*" end
  aligns = seperateLines(ali01, seperator_line_requirement)

  def alignParts(align)
    tt2 = align[0].split(' ')
    tt2[0][0, 1] = '' if tt2[0][0, 1] == '>'
    n1, n2 = 0, 0
    align_hash = Hash.new
    for xa in 0..align[1].size
      e1 = align[1][xa, 1]
      e2 = align[3][xa, 1]
      n1 += 1 if e1 != '-'
      n2 += 1 if e2 != '-'
      if e1 != '-' and e2 != '-'
        key = "#{tt2[0]}:#{n2}"
        align_hash[key] = n1
      end
    end

    return align_hash
  end

  all_aligns = aligns.collect { |align| alignParts(align) }.reduce(&:merge)

  ################################################################################
  # deal with paramsff
  parameters = File.open(paramsff).readlines.each do |line| line.chomp! end
  pmf_dist = Hash.new  # hash the pmf distance between atom pairs
  parameters.each do |line|
    if line.length > 4 and line[0, 3] == 'PMF'
      tt3 = line.split(' ')
      pmf_dist["#{tt3[1]}:#{tt3[2]}"] = tt3[3]
    end
  end


  ################################################################################
  # deal with pockets.dat
  pkt_lines = File.open(fpkt1).readlines.each do |line| line.chomp! end

  seperator_line_requirement = Proc.new do |line| line == 'TER' end

  pkts = seperateLines(pkt_lines, seperator_line_requirement)

  # only select pocket # equals fnum1 for one .pockets.dat file
  eligible_pkts = pkts.select { |pkt| pkt[0][7, 4].to_i == fnum1 }

  my_pkt = []
  if eligible_pkts.size > 0
    my_pkt  = eligible_pkts[0]
  else
    raise "NO eligible pockets found!!"
  end
  
  bind_residues = my_pkt.grep(/RESIDUE /).grep(/\*/)  # * stands for binding residues
  res_hash = Hash.new
  bind_residues.each do |res|
    res_num = res[8, 5].to_i
    res_hash[res_num] = res[39, 8].to_f  # last column
  end

  ################################################################################
  # deal with templates.pdb

  # mapping from three letters to one letter
  cod01 = { 'ALA' => 'A',
            'CYS' => 'C',
            'ASP' => 'D',
            'GLU' => 'E',
            'PHE' => 'F',
            'GLY' => 'G',
            'HIS' => 'H',
            'ILE' => 'I',
            'LYS' => 'K',
            'LEU' => 'L',
            'MET' => 'M',
            'ASN' => 'N',
            'PRO' => 'P',
            'GLN' => 'Q',
            'ARG' => 'R',
            'SER' => 'S',
            'THR' => 'T',
            'VAL' => 'V',
            'TRP' => 'W',
            'TYR' => 'Y' }


  # split one protein into its constructing amino acids
  def split2AminoAcids(atom_lines)
    atom_num_strs = atom_lines.collect { |line| line[22, 4] }.uniq
    amino_acids = []
    for str in atom_num_strs
      atoms = atom_lines.select { |line| line[22, 4] == str }
      amino_acids.push(atoms)
    end

    amino_acids
  end

  define_method :getPeptidePlaneCenter do |amino_acid|
    acid_name = amino_acid[0][17, 3]
    if acid_name == 'GLY'
      atom_num = amino_acid[0][22, 4].to_i.to_s
      atom_line = amino_acid.grep(/ CA /)[0]
      x = atom_line[30, 8].to_f.to_s
      y = atom_line[38, 8].to_f.to_s
      z = atom_line[46, 8].to_f.to_s
      return "G:#{atom_num}:#{x}:#{y}:#{z}&"
    elsif ['ALA', 'SER', 'THR', 'VAL', 'LEU', 'ILE', 'ASN', 'ASP', 'PRO', 'CYS'].include? acid_name
      atom_num = amino_acid[0][22, 4].to_i.to_s
      plane_atoms = amino_acid.reject { |atom| / N  | CA | C  | O  / =~ atom }
      if plane_atoms.size > 0
        x_coords, y_coords, z_coords = [], [], []
        plane_atoms.each do |atom_line|
          x_coords.push(atom_line[30, 8].to_f)
          y_coords.push(atom_line[38, 8].to_f)
          z_coords.push(atom_line[46, 8].to_f)
        end
        x_center = (x_coords.reduce(:+) / x_coords.size).to_s
        y_center = (y_coords.reduce(:+) / y_coords.size).to_s
        z_center = (z_coords.reduce(:+) / z_coords.size).to_s
        return "#{cod01[acid_name]}:#{atom_num}:#{x_center}:#{y_center}:#{z_center}&"
      else
        return ""
      end
    elsif ['ARG', 'LYS', 'GLU', 'GLN', 'HIS', 'MET', 'PHE', 'TYR', 'TRP'].include? acid_name
      atom_num = amino_acid[0][22, 4].to_i.to_s
      plane_atoms = amino_acid.reject { |atom|  / N  | CA | C  | O  | CB | CG / =~ atom }
      if plane_atoms.size > 0
        x_coords, y_coords, z_coords = [], [], []
        plane_atoms.each do |atom_line|
          x_coords.push(atom_line[30, 8].to_f)
          y_coords.push(atom_line[38, 8].to_f)
          z_coords.push(atom_line[46, 8].to_f)
        end
        x_center = (x_coords.reduce(:+) / x_coords.size).to_s
        y_center = (y_coords.reduce(:+) / y_coords.size).to_s
        z_center = (z_coords.reduce(:+) / z_coords.size).to_s
        return "#{cod01[acid_name]}:#{atom_num}:#{x_center}:#{y_center}:#{z_center}&"
      else
        return ""
      end
    end
  end

  pdb_lines = File.open(ftpl1, 'r').readlines.each do |line| line.chomp! end
  templates = Hash.new

  # seperate ecch protein by "TER"
  seperator_line_requirement = Proc.new do |line| line == "TER" end
  proteins = seperateLines(pdb_lines, seperator_line_requirement)

  templates = Hash.new
  proteins.each do |protein|
    protein_id = protein.grep(/REMARK   PDB-ID:/)[0].split(/\ /).pop
    effective_pts = ""
    atoms = protein.grep(/ATOM /)
    amino_acids = split2AminoAcids(atoms)
    amino_acids.each do |aa|
      effective_pts += getPeptidePlaneCenter(aa)
    end
    effective_pts[-1] = '' if effective_pts[-1] == '&'
    templates[protein_id] = effective_pts 
  end

  pmf_cnt = []
  eligible_ligs.each do |molid, lines|
    if templates.has_key? molid[0, 5]
      sdf_tmp = Tempfile.new("kde_sdf")
      path = sdf_tmp.path
      sdf_tmp.close
      File.open(path, 'w') do |file| file.write(lines.join("\n")) end

      # convert from sdf to mol2
      cmd = "#{babel} -isdf #{path} -omol2 -"
      stdout_str, stderr_str, status = Open3.capture3(cmd)
      unless status.success?
        STDERR.puts "Error running #{cmd}\n"
      end

      # grep atom lines
      mol2 = stdout_str.split("\n")
      lig_atoms = getLigAtomsFromMol2(mol2)
      sdf_tmp.unlink

      atom_strs = []
      lig_atoms.each do | atom | 
        atom_name, x, y, z = atom
        atom_strs.push("#{atom_name}:#{x}:#{y}:#{z}")
      end

      alignments = templates[molid[0, 5]].split(/\&/)
      alignments.each do |align|
        seqs = align.split(/\:/)
        pmf_key = "#{molid[0, 5]}:#{seqs[1]}"
        if all_aligns.has_key? pmf_key
          atom_strs.each do |atom|
            n0, x0, y0, z0 = atom.split(/:/)[0..3]
            x1, y1, z1 = seqs[2..4]
            dis = Math.sqrt((x0.to_f - x1.to_f)**2 + (y0.to_f - y1.to_f)**2 + (z0.to_f - z1.to_f)**2)
            atom_name = 'CA'

            atom_name = seqs[0] + '1' if seqs[0] != 'G'
            if pmf_dist.has_key? "#{atom_name}:#{n0}"
              pmf_val = pmf_dist["#{atom_name}:#{n0}"]
              pmf_cnt.push(all_aligns[pmf_key].to_s + ':' + n0)  if dis <= pmf_val.to_f and res_hash.has_key? all_aligns[pmf_key]
            end
          end
        end
      end
    end
  end

  num_pmf_cnt = pmf_cnt.size
  pmf_dist = Hash.new(0)
  pmf_cnt.each do | pmf |
    pmf_dist[pmf] += 1
  end

  kde_hash = Hash.new(0)
  kde_pts.each do |kde|
    n1 = kde.split(' ')[1]
    kde_hash[n1] += 1
  end

  atm01 = %w{ Br C.1 C.2 C.3 C.ar C.cat Cl F I N.1 N.2 N.3 N.4 N.am N.ar N.pl3 O.2 O.3 O.co2 P.3 S.2 S.3 S.O S.O2 };

  psp = []
  res_hash.keys.sort.each do |key|
    atm01.each do |atm|
      pmf1 = 0.0
      pmf1 = pmf_dist["#{key}:#{atm}"].to_f if pmf_dist.has_key? "#{key}:#{atm}"
      pmf2 = 0.0
      pmf2 = kde_hash[atm].to_f if kde_hash.has_key? atm
      pmf3 = 1.0 / bind_residues.size
      pmf4 = pmf2 / kde_pts.size
      pmf5 = num_pmf_cnt * pmf3 * pmf4

      if pmf1.to_f != 0.0 and pmf5.to_f != 0.0
        pmf6 = -1.0 * Math.log(pmf1 / pmf5)
        psp.push("PSP #{key} #{atm} #{pmf6}")
      end
    end
  end

  return psp
end


def prepareMCS(fsdf1, eligible_ligs, fkey1, pkcombu='')
  # tmp file for input sdf
  sdf_tmp = Tempfile.new("sdf")
  sdf_path = sdf_tmp.path
  sdf_tmp.close

  # tmp file eligible ligs
  eli_tmp = Tempfile.new("eli")
  eli_path = eli_tmp.path
  eli_tmp.close

  # tmp file for pkcombu output
  pk_tmp = Tempfile.new("pk")
  pk_path = pk_tmp.path
  pk_tmp.close

  # read sdf and seperate the ligands by $$$$
  sdf_lines = File.open(fsdf1).readlines.each { |line| line.chomp! }
  seperator_line_requirement = Proc.new do |line| line == "$$$$" end
  sdf_ligs = seperateLines(sdf_lines, seperator_line_requirement)

  # get ligand molid
  define_method :getMolid do |lig_lines|
    molid = ''
    lig_lines.each_with_index do |line, idx|
      molid = lig_lines[idx+1] if line =~ /#{fkey1}/
    end
    return molid
  end

  # compare the ligand with eligible templates
  mcs = []
  sdf_ligs.each do |lig_lines|
    molid = getMolid(lig_lines)
    if molid.size > 2
      puts molid + ' -> < MCS >'
      # write to tmp file for comparison
      File.open(sdf_path, 'w') do |file| file.write(lig_lines.join("\n")) end

      eligible_ligs.each do |eli_lig, eli_lines|
        puts "\t\t#{eli_lig}\n"
        mcs_line = ''

        File.open(eli_path, 'w') do |file| file.write(eli_lines.join("\n")) end

        # run pkcombu
        cmd = "#{pkcombu} -A #{sdf_path} -B #{eli_path} -oam #{pk_path} -fA S -fB S"
        stdout_str, stderr_str, status = Open3.capture3(cmd)
        cal1 = stdout_str.split("\n") if stdout_str.size > 0
        cal2 = File.open(pk_path).readlines.each { |line| line.chomp! }

        tcc = 0.0
        num_pairs = 0
        mapping = Hash.new
        cal2.each do |line|
          if line =~ /\#/ and line =~ /tanimoto/
            tcc = line.split[-1]
          elsif
            line =~ /^[0-9]/
            atomA = line.split[0]
            atomB = line.split[1]
            num_pairs += 1
            mapping[atomA] = atomB
          end
        end

        if num_pairs > 0
          eli_molid = getMolid(eli_lines)
          mcs_line = "MCS #{molid} #{eli_molid} #{tcc} #{num_pairs}"

          mapping.keys.sort.each do |key|
            val = mapping[key]
            x = eli_lines[val.to_i + 3][0, 10].to_f
            y = eli_lines[val.to_i + 3][10, 10].to_f
            z = eli_lines[val.to_i + 3][20, 10].to_f

            mcs_line += " #{key} #{x} #{y} #{z}"
          end
        end
        mcs.push(mcs_line)
      end
    end
  end

  sdf_tmp.unlink
  eli_tmp.unlink
  pk_tmp.unlink

  return mcs
end

