import numpy as np
import matplotlib.pyplot as plt

class Dane:
    czas_calkowity_TC = 2  # czas całkowity
    ilosc_bitow = 42
    docelowa_czestotliowsc_W = 2  # docelowa czestotliwosc
    A = 3
    h = 0

    czas_trwania_bitu_TB = czas_calkowity_TC / ilosc_bitow  # czas trwania pojedynczego bitu
    czestotliowsc_fali_nosnej_FN = docelowa_czestotliowsc_W * (1 / czas_trwania_bitu_TB)
    czestotliwosc_probkowania_FS = 16 * ilosc_bitow # czestotliwosc probkowania
    czas_t = np.arange(0, czas_calkowity_TC, 1 / czestotliwosc_probkowania_FS)
    Fn1 = (docelowa_czestotliowsc_W + 1) / czas_trwania_bitu_TB
    Fn2 = (docelowa_czestotliowsc_W + 2) / czas_trwania_bitu_TB
    probki_na_bit = len(czas_t) / ilosc_bitow

class Kluczowanie:
    dane = Dane()
    def ASCI_na_bity(self,ciag_znakow):
        result = ""
        for znak in ciag_znakow:
            ascii = ord(znak)
            if ascii >=32 and ascii <= 127:
                format_binarny = format(ascii,"07b")
                result += format_binarny
        return result


    probki_na_bit = int(len(dane.czas_t) / dane.ilosc_bitow)
    def kluczowanie_ASK(self,A1,A2,b):
        Za = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i/self.probki_na_bit)] == "0":
                Za.append(A1 * np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
            else:
                Za.append(A2 * np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
        return Za

    def kluczowanie_PSK(self,b):
        Zp = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i / self.probki_na_bit)] == "0":
                Zp.append(np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
            else:
                Zp.append(np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i] + np.pi))
        return Zp

    def kluczowanie_FSK(self,b):
        Zf = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i / self.probki_na_bit)] == "0":
                Zf.append(np.sin(2 * np.pi * self.dane.Fn1 * self.dane.czas_t[i]))
            else:
                Zf.append(np.sin(2 * np.pi * self.dane.Fn2 * self.dane.czas_t[i]))
        return Zf



    def oblicz_widmo(self,Zx):
        widmo = []
        widmo = np.abs(np.fft.fft(Zx))
        widmo = 10 * np.log10(widmo)
        return widmo

    def wyswietl_i_zapisz_widmo(self,x):
        nazwa = next((nazwa for nazwa, wartosc in globals().items() if wartosc is x), None)
        plt.title(f"{str(nazwa)}")
        plt.plot(x)
        plt.savefig(f"{str(nazwa)}.png")
        plt.show()

    def wyswietl_sygnal(self,z,nazwa, x = 1):
        plt.title(f"{nazwa}")
        if x == 1:
            plt.plot(z)
            plt.savefig(f"{nazwa}.png")
        else:
            plt.plot(z[:x])
            plt.savefig(f"{nazwa}.png")
        plt.show()

    def oszacuj_szerokosc_pasma(self,db,widmo):
        maxWidmo = np.max(widmo)
        prog = maxWidmo - db
        index = np.where(widmo > prog)[0]
        niska  = index[0]
        wysoka = index[-1]
        szerokosc = wysoka - niska
        print(szerokosc)


class Demodulacja:
    dane = Dane()
    def pomnoz(self,XSK,name,Fnx = True):
        if name == "ASK" or name == "PSK":
            xt = self.dane.A * np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t)
        else:
            if Fnx:
                xt = self.dane.A * np.sin(2 * np.pi * self.dane.Fn1 * self.dane.czas_t)
            else:
                xt = self.dane.A * np.sin(2 * np.pi * self.dane.Fn2 * self.dane.czas_t)
        return XSK * xt

    def podziel_sygnal_i_zsumuj(self,XSK):
        pt = []
        tmp = 0
        for i in range(0, len(XSK), int(self.dane.probki_na_bit)):
            for j in range(int(self.dane.probki_na_bit)):
                tmp += XSK[i + j]
            pt.append(tmp)
            tmp = 0
        self.ustaw_h(pt)
        return pt

    def  dzialanie_na_p1_p2_FSK(self,p1,p2):
        wynik = []
        for a, b in zip(p2, p1):
            wynik.append(a - b)
        return wynik

    def ustaw_h(self,pt):
        self.dane.h = np.average(pt)

    def zwroc_bitowo(self,XSK,x):
        ct = []
        if x == "ASK":
            for i in XSK:
                if i > self.dane.h:
                    ct.append(1)
                else:
                    ct.append(0)
        elif x == "PSK":
            for i in XSK:
                if i < 0:
                    ct.append(1)
                else:
                    ct.append(0)
        else:
            for i in XSK:
                if i > 0:
                    ct.append(1)
                else:
                    ct.append(0)

        return ct

    def zamien_na_bity(self,ct):
        result = ""
        for i in ct:
            result += str(i)
        return result






kluczowanie = Kluczowanie()
demodulacja = Demodulacja()
ASK = "ASK"
PSK = "PSK"
FSK = "FSK"
ciag = kluczowanie.ASCI_na_bity("Robert")

print("DEMODULACJA ASK:")
print(f"Przed: {ciag}")
sygnal_ASK = kluczowanie.kluczowanie_ASK(1, 2, ciag)
xt = demodulacja.pomnoz(sygnal_ASK,ASK)
pt = demodulacja.podziel_sygnal_i_zsumuj(xt)
ret = demodulacja.zwroc_bitowo(pt,ASK)
print(f"PO:    {demodulacja.zamien_na_bity(ret)}\n\n")

kluczowanie.wyswietl_sygnal(sygnal_ASK, "ask_z")
kluczowanie.wyswietl_sygnal(xt,         "ask_x")
kluczowanie.wyswietl_sygnal(pt,         "ask_p")
kluczowanie.wyswietl_sygnal(ret,        "ask_c")


print("DEMODULACJA PSK:")
print(f"Przed: {ciag}")
sygnal_PSK = kluczowanie.kluczowanie_PSK(ciag)
xt = demodulacja.pomnoz(sygnal_PSK,PSK)
pt = demodulacja.podziel_sygnal_i_zsumuj(xt)
ret = demodulacja.zwroc_bitowo(pt,PSK)
print(f"PO:    {demodulacja.zamien_na_bity(ret)}\n\n")

kluczowanie.wyswietl_sygnal(sygnal_PSK, "psk_z")
kluczowanie.wyswietl_sygnal(xt,         "psk_x")
kluczowanie.wyswietl_sygnal(pt,         "psk_p")
kluczowanie.wyswietl_sygnal(ret,        "psk_c")

# FNx = true  -> FN1
# FNx = false -> FN2

print("DEMODULACJA FSK:")
print(f"Przed: {ciag}")
sygnal_FSK = kluczowanie.kluczowanie_FSK(ciag)
xt_FN1 = demodulacja.pomnoz(sygnal_FSK,FSK)
xt_FN2 = demodulacja.pomnoz(sygnal_FSK,FSK,False)

pt_FN1 = demodulacja.podziel_sygnal_i_zsumuj(xt_FN1)
pt_FN2 = demodulacja.podziel_sygnal_i_zsumuj(xt_FN2)

pt =  demodulacja.dzialanie_na_p1_p2_FSK(pt_FN1,pt_FN2)
ret = demodulacja.zwroc_bitowo(pt,FSK)

kluczowanie.wyswietl_sygnal(sygnal_FSK, "fsk_z")
kluczowanie.wyswietl_sygnal(xt_FN1, "fsk_x1")
kluczowanie.wyswietl_sygnal(xt_FN2, "fsk_x2")
kluczowanie.wyswietl_sygnal(pt_FN1, "fsk_p1")
kluczowanie.wyswietl_sygnal(pt_FN2, "fsk_p2")
kluczowanie.wyswietl_sygnal(ret,     "fsk_c")

print(f"PO:    {demodulacja.zamien_na_bity(ret)}\n\n")





