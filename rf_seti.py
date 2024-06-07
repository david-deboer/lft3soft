from numpy import pi

class FMStation:
  def __init__(self, freq, mod, Tx, distance):
    self.freq = freq
    self.mod = mod
    self.power = Tx / (4.0 * pi * distance**2)

krfi = FMStation(freq=0.7e6, mod=100.0e3, Tx=1000, distance=50000)
wtfr = FMStation(freq=1.1e6, mod=110.0e3, Tx=2000, distance=35000)


class Location:
  def __init__(self, Tsys, D, delay, fm):
    self.Tsys = Tsys
    self.A = pi * (D / 2.0)**2
    self.site_delay = delay
    self.fm = fm
  
site1 = Location(Tsys=30.0, D=40.0, delay=0.0, fm=krfi)
site2 = Location(Tsys=30.0, D=40.0, delay=0.0, fm=wtfr)