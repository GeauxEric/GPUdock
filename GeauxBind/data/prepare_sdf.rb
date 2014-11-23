def addMolid(ifn, ofn)
  lines = File.open(ifn).readlines
  molid = lines[0].chomp
  to_add = "\n>  <MOLID> (#{molid})\n#{molid}\n\n"
  lines.insert(-2, to_add)
  File.open(ofn, 'w') {|file| file.write(lines.join)}
end

  
if ARGV[0] == __FILE__
  ifn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158.sdf1"
  ofn = "/home/jaydy/Workspace/GitHub/GPUdock/GeauxBind/data/ZINC00002158.sdf2"
  addMolid(ifn, ofn)
