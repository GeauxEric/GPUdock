require 'test/unit'
require 'open3'


class PrepareSdfTest < Test::Unit::TestCase

  def test_model_ens
    work_dir = "../data"
    prt_code = '1b9vA'

    py = "../src/run_prt_ens.py"
    cmd = "cd #{work_dir} && python #{py} -d . -p #{prt_code}"

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
