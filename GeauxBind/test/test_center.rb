require 'test/unit'
require_relative '../src/ff'

class CenterTest < Test::Unit::TestCase
  def test_a_getPocketCenterCoords
    fpkt1 = "../data/1b9vA.pockets.dat"
    fnum1 = 1
    center = getPocketCenterCoords(fpkt1, fnum1)
    assert_equal('31.860', center.split[1])
    assert_equal('-9.590', center.split[2])
    assert_equal('64.299', center.split[3])
  end

end
