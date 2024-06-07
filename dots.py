import matplotlib.pyplot as plt
import numpy as np

astrobotics = 3.5
firefly = 103 * 2.54 / 100
im = 4.0

r = astrobotics / 2.0
spacing = 0.5
extra = {'lo': [0.0, 1.5], 'fm': [1.5, 0.0], 'hi': [-1.5, 0.0], 'comm': [0.0, -1.5]}

x = []
y = []

for th in np.linspace(0.0, 2.0 * np.pi, 360 * 2):
    x.append(r * np.cos(th))
    y.append(r * np.sin(th))
plt.plot(x, y)

dox = []
doy = []

remove = [0.25, 1.25]

num = 0
for i in np.arange(-r, r + spacing, spacing):
    for j in np.arange(-r, r + spacing, spacing):
        this_r = np.sqrt(i*i + j*j)
        if abs(i) in remove and abs(j) in remove and this_r > 1.0:
            continue
        if this_r <= r:
            num += 1
            dox.append(i)
            doy.append(j)
plt.plot(dox, doy, 'o', label='array')

for nm, po in extra.items():
    plt.plot(po[0], po[1], 'o', label=nm)

print(f"Total number = {num}")
plt.axis('image')
plt.legend(bbox_to_anchor=(1.3, 1.01))