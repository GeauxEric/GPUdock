require 'pry'
require 'tempfile'
require 'open3'

# use the requirement_meet method to seperate the text into several sections
def splitLines(lines, requirement_meet)
  splitted = []
  section = []
  lines.each do |line| 
    section.push(line)
    if requirement_meet.call(line)
      splitted.push(section)
      section = []
    end
  end

  splitted
end

# load a .ligands.sdf file, get the ligands with required pocket number
def getEligibleLigs(flig1, fnum1)

  def getPocketNum(lig)
    pocket = 0
    lig.each_with_index do |line, idx| 
      if line =~ /EFINDSITE_POCKET/
        pocket = lig[idx + 1]
      end
    end

    pocket.to_i
  end

  def getMolid(lig)
    molid = ''
    lig.each_with_index do |line, idx| 
      if line =~ /MOLID/
        molid = lig[idx + 1]
      end
    end
    molid
  end

  def requiredPocketNum?(lig, fnum)
    getPocketNum(lig) == fnum ? true : false
  end

  lig01 = File.readlines(flig1).each{ |line| line.chomp! }

  seperator_line_requirement = Proc.new do |line| 
    line == "$$$$"
  end

  ligs = splitLines(lig01, seperator_line_requirement)

  eligible_ligs = Hash.new
  ligs.each do |lig| 
    if requiredPocketNum?(lig, fnum1)
      molid = getMolid(lig)
      eligible_ligs[molid] = lig
    end
  end
  
  return eligible_ligs
end

# convert the sdf to mol2 format
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

  atom_line_is_found = false

  kde_atom_num = 0
  kde_pts = []
  mol2.each do |line|
    atom_line_is_found = false if line == '@<TRIPOS>BOND'

    if atom_line_is_found
      kde_atom_num += 1
      x1 = line[16, 10].to_f.to_s
      y1 = line[26, 10].to_f.to_s
      z1 = line[36, 10].to_f.to_s
      n1 = line[47, 7].delete(' ')

      kde_pts.push('KDE ' + [n1, x1, y1, z1].join(' ') + "\n")
    end
    atom_line_is_found = true if ( line == '@<TRIPOS>ATOM' );
  end

  return [kde_pts, kde_atom_num]
end

def preparePSP(fali1, paramsff, fpkt1, ftpl1, eligible_ligs, fnum1, kde_pts, kde_atom_num, babel='')
  # deal with .alignments.dat
  ali01 = File.open(fali1, 'r').readlines.each do |line| line.chomp! end
  seperator_line_requirement = Proc.new do |line| 
    line == "*"
  end
  aligns = splitLines(ali01, seperator_line_requirement)


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
        key = "#{tt2[0]}:#{n2}"  # what is this ???
        align_hash[key] = n1
      end
    end

    return align_hash
  end

  all_aligns = aligns.collect { |align| alignParts(align) }.reduce(&:merge)

  # deal with paramsff
  dat01 = File.open(paramsff).readlines.each do |line| line.chomp! end
  dat02 = Hash.new  # hash the pmf distance between atom pairs
  dat01.each do |line| 
    if line.length > 4 and line[0, 3] == 'PMF'
      tt3 = line.split(' ')
      dat02["#{tt3[1]}:#{tt3[2]}"] = tt3[3]
    end
  end


  # deal with pockets.dat
  # store the ligands with pockets number == to fnum1
  pkt01 = File.open(fpkt1).readlines.each do |line| line.chomp! end 
  pkt02, pkt03 = [], []

  pkt01.each do |line01| 
    pkt03.push(line01)
    if line01 == 'TER' and pkt03[0][7, 4].to_i == fnum1  # check POCKET
      pkt02 = []
      pkt03.each do |line03| 
        pkt02.push(line03)
      end
      pkt03 = []
    end
  end

  residues = pkt02.grep(/RESIDUE /).grep(/\*/)  # what doe star residues stand for?
  res_hash = Hash.new
  residues.each do |res| 
    tt4 = res[8, 5].to_i
    res_hash[tt4] = res[39, 8].to_f  # last column
  end

  # deal with templates.pdb
  pdb01 = File.open(ftpl1, 'r').readlines.each do |line| line.chomp! end
  pdb02 = Hash.new
  pdb03 = []

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

  def split2AminoAcids(atom_lines)
    atom_num_strs = atom_lines.collect { |line| line[22, 4] }.uniq
    amino_acids = []
    for str in atom_num_strs
      atoms = atom_lines.select { |line| line[22, 4] == str }
      amino_acids.push(atoms)
    end

    amino_acids
  end

  define_method :getNoneBackboneCenter do |amino_acid| 
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
      none_backbones = amino_acid.reject { |atom| / N  | CA | C  | O  / =~ atom }
      if none_backbones.size > 0
        x_coords, y_coords, z_coords = [], [], []
        none_backbones.each do |atom_line| 
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
      none_backbones = amino_acid.reject { |atom|  / N  | CA | C  | O  | CB | CG / =~ atom }
      if none_backbones.size > 0
        x_coords, y_coords, z_coords = [], [], []
        none_backbones.each do |atom_line| 
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


    seperator_line_requirement = Proc.new do |line| 
      line == "TER"
    end

    proteins = splitLines(pdb01, seperator_line_requirement)

    pdb02 = Hash.new
    proteins.each do |protein| 
      protein_id = protein.grep(/REMARK   PDB-ID:/)[0].split(/\ /).pop
      pdb04= ""
      atoms = protein.grep(/ATOM /)
      amino_acids = split2AminoAcids(atoms)
      amino_acids.each do |aa| 
        pdb04 += getNoneBackboneCenter(aa)
      end
      pdb04[-1] = '' if pdb04[-1] == '&'
      pdb02[protein_id] = pdb04
    end


    pmf_cnt = []
    eligible_ligs.each do |molid, lines| 
      if pdb02.has_key? molid[0, 5]
        sdf_tmp = Tempfile.new("kde_sdf")
        path = sdf_tmp.path
        sdf_tmp.close
        File.open(path, 'w') do |file|
          file.write(lines.join("\n"))
        end
        
        # convert from sdf to mol2
        cmd = "#{babel} -isdf #{path} -omol2 -"
        stdout_str, stderr_str, status = Open3.capture3(cmd)
        unless status.success?
          STDERR.puts "Error running #{cmd}\n"
        end


        # grep atom lines
        mol2 = stdout_str.split("\n")
        mol2_atoms = []
        atom_line_is_found = false
        mol2.each do |line| 
          atom_line_is_found = false if line == '@<TRIPOS>BOND'
          if atom_line_is_found
            x = line[16, 10].to_f.to_s
            y = line[26, 10].to_f.to_s
            z = line[36, 10].to_f.to_s
            atom_name = line[47, 7].delete(' ')
            mol2_atoms.push("#{atom_name}:#{x}:#{y}:#{z}")
          end
          atom_line_is_found = true if line == '@<TRIPOS>ATOM'
        end
        sdf_tmp.unlink

        alignments = pdb02[molid[0, 5]].split(/\&/)
        alignments.each do |align| 
          seqs = align.split(/\:/)
          pmf_key = "#{molid[0, 5]}:#{seqs[1]}"
          if all_aligns.has_key? pmf_key
            mol2_atoms.each do |atom| 
              n0, x0, y0, z0 = atom.split(/:/)[0..3]
              x1, y1, z1 = seqs[2..4]
              dis = Math.sqrt((x0.to_f - x1.to_f)**2 + (y0.to_f - y1.to_f)**2 + (z0.to_f - z1.to_f)**2)
              atom_name = 'CA'

              atom_name = seqs[0] + '1' if seqs[0] != 'G'
              if dat02.has_key? "#{atom_name}:#{n0}"
                pmf_val = dat02["#{atom_name}:#{n0}"]
                pmf_cnt.push(all_aligns[pmf_key].to_s + ':' + n0)  if dis <= pmf_val.to_f and res_hash.has_key? all_aligns[pmf_key]
              end
            end
          end
        end
      end
    end

    num_pmf_cnt = pmf_cnt.size
    pmf_hash = Hash.new(0)
    pmf_cnt.each do | pmf |
      pmf_hash[pmf] += 1
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
        pmf1 = pmf_hash["#{key}:#{atm}"].to_f if pmf_hash.has_key? "#{key}:#{atm}"
        pmf2 = 0.0
        pmf2 = kde_hash[atm].to_f if kde_hash.has_key? atm
        pmf3 = 1.0 / residues.size
        pmf4 = pmf2 / kde_atom_num
        pmf5 = num_pmf_cnt * pmf3 * pmf4

        if pmf1.to_f != 0.0 and pmf5.to_f != 0.0
          pmf6 = -1.0 * Math.log(pmf1 / pmf5)
          psp.push("PSP #{key} #{atm} #{pmf6}")
        end
      end
    end

    return psp
end

