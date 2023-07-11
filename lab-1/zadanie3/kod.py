import numpy as np
import matplotlib.pyplot as plot

tc = 10       # czas calkowity w sekundach
fs = 8000     # czestotliwosc/okres próbkowania HZ ile probek na sekunde
N  = tc * fs  # liczba probek przypadajacych na caly sygnal
fi = 0        # faza sygnalu
ts = 1/fs     # okres probkowania
f  = 5

t = np.arange(0, tc, ts) #Wygenerowane probki
u = []

# Tabela 3 funkcja 3
for time in t:
    if time < 1.2 and time >= 0:
        u.append((-(time * time) + 0.5) * np.sin(30 * np.pi * time) * np.log2(time * time + 1))

    elif time < 2 and time >=1.2:
        u.append((1/time) * 0.8 * np.sin(24 * np.pi * time) - 0.1 * time)
    elif time < 2.4 and time >= 2:
        u.append(np.abs(np.sin(2 * np.pi * time * time)) ** 0.8)
    else:
        u.append(0.23 * np.sin(20 * np.pi * time) * np.sin(12 * np.pi * time))

plot.xlabel("t")
plot.ylabel("x(t)")
plot.title("Franczak uZadanie 3, funkcja nr 2")
plot.plot(t, u)
plot.show()

