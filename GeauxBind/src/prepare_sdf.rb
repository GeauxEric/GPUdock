# require 'pry'
require 'tempfile'


# concat the content in ifns into ofn
def concatFiles(ifns, ofn)
  ifns.each do
    |ifn| raise "#{ifn} does not exist!" unless File.exist?(ifn)
  end
  to_concat = ifns.join(' ')
  cmd = "cat #{to_concat} > #{ofn}"
  concat_well = system cmd
  raise "failure to concat #{to_concat}" unless concat_well
end


# remove the hydrogen atom by using babel
def removeH(ifn, ofn, babel='')
  if File::exists?(babel) and File::exists?(ifn)
    babel_succeed = system "#{babel} -isdf #{ifn} -osdf #{ofn} -d"
    raise "failure in running babel" unless babel_succeed
  else
    raise "certain files do not exist!\n"
  end

end


def addMolid(lines, key="$$$$")
  molid = lines[0].chomp
  raise "invalid sdf file!" if molid.length < 2

  added_lines = []
  lines.each_with_index {
    |line, idx|
    if line.include? key
      to_add = "\n>  <MOLID> (#{molid})\n#{molid}\n\n$$$$\n"
      added_lines.push(to_add)
      # if not end of file, next line supposed to contain another molecure id
      molid = lines[idx + 1].chomp if idx + 1 < lines.size
    else
      added_lines.push(line)
    end
  }
  return added_lines
end


# add molecule ID to the sdf file
def addMolid2File(ifn, ofn)
  if File::exists?(ifn)
    lines = File.open(ifn).readlines
    added_lines = addMolid(lines, key="$$$$")
    File.open(ofn, 'w') {|file| file.write(added_lines.join)}
  else
    raise "#{ifn} does not exist\n"
  end

end
