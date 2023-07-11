import matplotlib.pyplot as plt
import numpy as np

#----------------
Tc = 1
fs = 1000
Ts = 1 / fs
N = Tc * fs
t = np.arange(0, Tc, Ts)
f = 5
y = np.cos(2 * np.pi * f * t) * np.cos(2.5 * t ** 0.2 * np.pi)
#----------------
dftI = 0
dft = []
widmo = []
widmoDec = []
skalaCzest = []
#----------------
for i in range(int(N / 2)):
    for n in range(N):
        dftI += y[n] * np.exp(-2j * np.pi * i * n / N)
    dft.append(dftI)
    widmo.append(np.sqrt(dft[i].real ** 2 + dft[i].imag ** 2))
    widmoDec.append(10 * np.log(widmo[i]))
    skalaCzest.append(i * fs / N)

plt.figure(1)
plt.xlabel("[Hz]")
plt.ylabel("[dec]")
plt.stem(skalaCzest, widmoDec)
plt.xlim(0, 10)
plt.savefig("dft.png")
plt.show()

plt.figure(2)
plt.xlabel("s]")
plt.ylabel("y")
plt.plot(t, y)
plt.savefig("y.png")
plt.show()
