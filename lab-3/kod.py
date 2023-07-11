import numpy as np
import matplotlib.pyplot as plt

#Czestotliwosci
fm = 2
fn = 50

#zmienne do czasu
Tc = 1
czestotliwoscProbkowania = 8000

#CZAS
t = np.arange(0,Tc,1/czestotliwoscProbkowania)
mt = np.sin(2 * np.pi * fm * t)

#Zadanie 1

# A) Modulacja amplitudy
def obliczKA(Ka):
    return (Ka * mt + 1) * np.cos(2 * np.pi * fn * t)

kaA = obliczKA(0.5)
kaB = obliczKA(6)
kaC = obliczKA(22)

# B) Modulacja fazy
def obliczKP(Kp):
    return np.cos(2 * np.pi * fn * t + Kp * mt)

kpA = obliczKP(0)
kpB = obliczKP(2.40)
kpC = obliczKP(3* np.pi)

# C) Modulacja czestotliwosci
def obliczKF(Kf):
    return np.cos(2 * np.pi * fn * t + (Kf/fn) * mt)

kfA = obliczKF(0)
kfB = obliczKF(2.40)
kfC = obliczKF(3* np.pi)

# !!!!!!!!!!!!TU PAN PODAL NA TABLICY TO SIE WSPOMAGAM :D
#///////// ZA
za_a = np.fft.fft(kaA)
za_a = np.abs(za_a)
za_a = 10 * np.log10(za_a)[0:60]
plt.plot(za_a)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("za_a.png")
plt.show()

za_b = np.fft.fft(kaB)
za_b = np.abs(za_b)
za_b = 10 * np.log10(za_b)[0:60]
plt.plot(za_b)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("za_b.png")
plt.show()

za_c = np.fft.fft(kaC)
za_c = np.abs(za_c)
za_c = 10 * np.log10(za_c)[0:60]
plt.plot(za_c)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("za_c.png")
plt.show()

# ///////////////////  ZP
zp_a = np.fft.fft(kpA)
zp_a = np.abs(zp_a)
zp_a = 10 * np.log10(zp_a)[0:100]
plt.plot(zp_a)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zp_a.png")
plt.show()

zp_b = np.fft.fft(kpB)
zp_b = np.abs(zp_b)
zp_b = 10 * np.log10(zp_b)[0:100]
plt.plot(zp_b)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zp_b.png")
plt.show()

zp_c = np.fft.fft(kpC)
zp_c = np.abs(zp_c)
zp_c = 10 * np.log10(zp_c)[0:100]
plt.plot(zp_c)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zp_c.png")
plt.show()

#/////////////////// ZF

zf_a = np.fft.fft(kfA)
zf_a = np.abs(zf_a)
zf_a = 10 * np.log10(zf_a)[0:100]
plt.plot(zf_a)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zf_a.png")
plt.show()

zf_b = np.fft.fft(kfB)
zf_b = np.abs(zf_b)
zf_b = 10 * np.log10(zf_b)[0:100]
plt.plot(zf_b)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zf_b.png")
plt.show()

zf_c = np.fft.fft(kfC)
zf_c = np.abs(zf_c)
zf_c = 10 * np.log10(zf_c)[0:100]
plt.plot(zf_c)
plt.xlabel('indeks widma')
plt.ylabel('amplituda widma')
plt.title('widmo f=1 kHz')
plt.savefig("zf_c.png")
plt.show()