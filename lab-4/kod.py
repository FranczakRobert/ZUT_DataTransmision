import numpy as np
import matplotlib.pyplot as plt

class Dane:
    czas_calkowity_TC = 1  # czas całkowity
    ilosc_bitow = 10
    docelowa_czestotliowsc_W = 2  # docelowa czestotliwosc

    czas_trwania_bitu_TB = czas_calkowity_TC / ilosc_bitow  # czas trwania pojedynczego bitu
    czestotliowsc_fali_nosnej_FN = docelowa_czestotliowsc_W * (1 / czas_trwania_bitu_TB)
    czestotliwosc_probkowania_FS = 16 * czestotliowsc_fali_nosnej_FN # czestotliwosc probkowania
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

    def ustaw_dane(self):
        dane = Dane()

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

    def wyswietl_sygnal(self, z, nazwa, x=1):
        plt.title(f"{nazwa}")
        if x == 1:
            plt.plot(self.dane.czas_t, z)
            # plt.savefig(f"{nazwa}.png")
        else:
            plt.plot(self.dane.czas_t[:x], z[:x])
            # plt.savefig(f"{nazwa}.png")
        plt.show()

    def oszacuj_szerokosc_pasma(self,db,widmo):
        maxWidmo = np.max(widmo)
        prog = maxWidmo - db
        index = np.where(widmo > prog)[0]
        niska  = index[0]
        wysoka = index[-1]
        szerokosc = wysoka - niska
        print(szerokosc)


x = Kluczowanie()

b = x.ASCI_na_bity("Robert Franczak")[:10]
Za = x.kluczowanie_ASK(1,2,b)
Zp = x.kluczowanie_PSK(b)
Zf = x.kluczowanie_FSK(b)

x.wyswietl_sygnal(Za, "za")
x.wyswietl_sygnal(Zp, "zp")
x.wyswietl_sygnal(Zf, "zf")

za_widmo = x.oblicz_widmo(Za)
zp_widmo = x.oblicz_widmo(Zp)
zf_widmo = x.oblicz_widmo(Zf)

x.wyswietl_i_zapisz_widmo(za_widmo)
x.wyswietl_i_zapisz_widmo(zp_widmo)
x.wyswietl_i_zapisz_widmo(zf_widmo)
print(" ")

db = 3
print(f"{db} db dla ZA,ZP,ZF:")
x.oszacuj_szerokosc_pasma(db, za_widmo)
x.oszacuj_szerokosc_pasma(db, zp_widmo)
x.oszacuj_szerokosc_pasma(db, zf_widmo)
print(" ")
db = 6
print(f"{db} db dla ZA,ZP,ZF:")
x.oszacuj_szerokosc_pasma(db, za_widmo)
x.oszacuj_szerokosc_pasma(db, zp_widmo)
x.oszacuj_szerokosc_pasma(db, zf_widmo)
print(" ")
db = 12
print(f"{db} db dla ZA,ZP,ZF:")
x.oszacuj_szerokosc_pasma(db, za_widmo)
x.oszacuj_szerokosc_pasma(db, zp_widmo)
x.oszacuj_szerokosc_pasma(db, zf_widmo)