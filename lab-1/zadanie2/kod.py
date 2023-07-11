import numpy as np
import matplotlib.pyplot as plot

tc = 10       # czas calkowity w sekundach
fs = 8000     # czestotliwosc/okres próbkowania HZ ile probek na sekunde
N  = tc * fs  # liczba probek przypadajacych na caly sygnal
fi = 0        # faza sygnalu
ts = 1/fs     # okres probkowania
f  = 5

#Od 0 do CZAS podzielone na kawalki
t = np.arange(0, tc, ts) #Wygenerowane probki

# FUNKCJA 1 Z TABELI 1
x = np.cos(2 * np.pi * f * t) * np.cos(2.5 * t ** 0.2 * np.pi)

# Zestaw 7
y = np.sin(np.pi * t) * np.sin(2 * x * np.pi * x)
plot.xlabel("Czas w sekundach")
plot.title("Franczak: zadanie2 - FUNKCJA: y")
plot.plot(t,y)
plot.show()

z = np.sqrt(np.abs(y)) - 3 * x
plot.xlabel("Czas w sekundach")
plot.title("Franczak: zadanie2 - FUNKCJA: z")
plot.plot(t,z)
plot.show()

v = x * y * y - z * np.cos(x)
plot.xlabel("Czas w sekundach")
plot.title("Franczak: zadanie2 - FUNKCJA: v")
plot.plot(t,v)
plot.show()