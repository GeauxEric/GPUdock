# require 'pry'

def removeH(ifn, ofn, babel='')
  # remove the hydrogen atom by using babel
  if File::exists?(babel)
    babel_succeed = system "#{babel} -isdf #{ifn} -osdf #{ofn} -d"
    raise "failure in running babel" unless babel_succeed
  else
    puts "please provide valide babel bin path"
  end

end



def addMolid(ifn, ofn)
  if File::exists?(ifn)
    lines = File.open(ifn).readlines
    molid = lines[0].chomp
    to_add = "\n>  <MOLID> (#{molid})\n#{molid}\n\n"
    lines.insert(-2, to_add)
    File.open(ofn, 'w') {|file| file.write(lines.join)}
  else
    puts "please provide a valid sdf file"
  end
  
end

