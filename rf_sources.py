from numpy import pi

class FMStation:
  def __init__(self, freq, mod, Tx, distance):
    self.freq = freq
    self.mod = mod
    self.power = Tx / (4.0 * pi * distance**2)

krfi = FMStation(freq=0.7e6, mod=100.0e3, Tx=1000, distance=50000)
kwtf = FMStation(freq=1.1e6, mod=110.0e3, Tx=2000, distance=35000)
