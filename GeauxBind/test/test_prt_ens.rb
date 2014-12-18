require 'test/unit'
require 'open3'


class PrepareSdfTest < Test::Unit::TestCase

  def test_model_ens
    work_dir = "."
    prt_code = '1b9vA'


    cmd = "python ../src/model_ens.py -d #{work_dir} -p #{prt_code}"
    Dir.chdir(work_dir)

    puts "\nRunning\t\t#{cmd}\n"
    stdout_str, stderr_str, status = Open3.capture3(cmd)

    # if status.success?
    #   puts "write to\t#{ofn}\n"
    # else
    #   STDERR.puts "Error running #{cmd}\n"
    #   exit 1
    # end
  end


end
