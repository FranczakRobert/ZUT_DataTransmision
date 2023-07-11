import numpy as np
import time
#---------
Tc = 1
fs = 1000
Ts = 1 / fs
N = Tc * fs #liczba próbek sygnalów w dziedzinie czasu i czestotliwosci.
t = np.arange(0, Tc, Ts)
f = 5
z = 0
y = np.cos(2 * np.pi * f * t) * np.cos(2.5 * t ** 0.2 * np.pi)
#---------

dft = []
start = time.time()
for i in range(N):
    for n in range(N):
        z += y[n] * np.exp(-2j * np.pi * i * n / N)
    dft.append(z)

end = time.time()
print("DFT TIME:") 
print(end - start)
