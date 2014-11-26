# require 'pry'

def removeH(ifn, ofn, babel='')
  # remove the hydrogen atom by using babel
  if File::exists?(babel) and File::exists?(ifn)
    babel_succeed = system "#{babel} -isdf #{ifn} -osdf #{ofn} -d"
    raise "failure in running babel" unless babel_succeed
  else
    raise "certain files do not exist!\n"
  end

end


def addMolid(lines)
  molid = lines[0].chomp
  to_add = "\n>  <MOLID> (#{molid})\n#{molid}\n\n"
  lines.insert(-2, to_add)
  return lines
end


def addMolid2File(ifn, ofn)
  if File::exists?(ifn)
    lines = File.open(ifn).readlines
    added_lines = addMolid(lines)
    File.open(ofn, 'w') {|file| file.write(added_lines.join)}
  else
    raise "#{ifn} does not exist\n"
  end
  
end

def addMultipleMolids(ifns, ofn)
  ifns.each {
    |ifn|
    if File::exists?(ifn)
      lines = File.open(ifn).readlines
      added_lines = addMolid(lines)
      File.open(ofn, 'a') {|file| file.write(added_lines.join)}
    else
      raise "#{ifn} does not exist"
    end
  }
end
