import numpy as np
import matplotlib.pyplot as plot

tc = 1        # czas calkowity w sekundach
fs = 22050     # czestotliwosc/okres próbkowania HZ ile probek na sekunde
N  = tc * fs  # liczba probek przypadajacych na caly sygnal
fi = 0        # faza sygnalu
ts = 1/fs     # okres probkowania

t = np.arange(0, tc, ts) #Wygenerowane probki

sum = []
# Tabela 4 funkcja 11
K = [1, 2, 3]
H = [2, 4, 16]
h = 0

for k in K:
    plot.xlabel("t")
    plot.ylabel("x(t)")
    plot.title("Zadanie 4, B" + str(h+1))
    print(H[h])
    sum.append((np.cos(12 * t * H[h] * H[h]) + np.cos(16 * t * H[h]))/(H[h] * H[h]))
    result = np.array(sum).sum(axis=0)
    print(result)
    h = h + 1
    plot.plot(t, result)
    plot.show()

    # Mam dosc...