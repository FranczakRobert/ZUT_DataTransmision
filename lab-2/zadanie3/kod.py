import matplotlib.pyplot as plt
import numpy as np
import time

#------------------------------
Tc = 1
fs = 1000
N = Tc * fs
Ts = 1 / fs
t = np.arange(0, Tc, Ts)
y = np.sin(2 * np.pi * 2 * t)
#------------------------------
dftI = 0
dft = []
widmo = []
widmoDec = []
skalaCzest = []
#------------------------------
start = time.time()
fft = np.fft.fft(y)
end = time.time()
print("Czas FFT:")
print(end - start)
#------------------------------
for i in range(int(N / 2)):
    widmo.append(np.sqrt(fft[i].real ** 2 + fft[i].imag ** 2))
    widmoDec.append(10 * np.log(widmo[i]))
    skalaCzest.append(i * fs / N)

plt.xlabel("[Hz]")
plt.ylabel("[dec]")
plt.stem(skalaCzest, widmoDec)
plt.xlim(0, 10)
plt.savefig("fft.png")
plt.show()
