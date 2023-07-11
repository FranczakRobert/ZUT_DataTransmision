import numpy as np
import matplotlib.pyplot as plot

tc = 10       # czas calkowity w sekundach
fs = 8000     # czestotliwosc/okres próbkowania HZ ile probek na sekunde
N  = tc * fs  # liczba probek przypadajacych na caly sygnal
fi = 0        # faza sygnalu
ts = 1/fs     # okres probkowania
f  = 5

# to wychodzi ze fs > 2 * fmax. Wiec fmax < 4000

#Od 0 do CZAS podzielone na kawalki
t = np.arange(0, tc, ts) #Wygenerowane probki

# FUNKCJA 1 Z TABELI 1
x = np.cos(2 * np.pi * f * t) * np.cos(2.5 * t ** 0.2 * np.pi)

plot.xlabel("Czas w sekundach")
plot.title("Franczak: zadanie1")
plot.plot(t,x)
plot.show()