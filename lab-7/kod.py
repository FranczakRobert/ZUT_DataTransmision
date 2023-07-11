import numpy as np
import matplotlib.pyplot as plt

class Dane:
    czas_calkowity_TC = 2  # czas całkowity
    ilosc_bitow = 196
    docelowa_czestotliowsc_W = 2  # docelowa czestotliwosc
    A = 3
    h = 0

    czas_trwania_bitu_TB = czas_calkowity_TC / ilosc_bitow  # czas trwania pojedynczego bitu
    czestotliowsc_fali_nosnej_FN = docelowa_czestotliowsc_W * (1 / czas_trwania_bitu_TB)
    czestotliwosc_probkowania_FS = 16 * ilosc_bitow # czestotliwosc probkowania
    czas_t = np.arange(czestotliwosc_probkowania_FS * czas_calkowity_TC) / czestotliwosc_probkowania_FS
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
    def zmien_strina_na_tablice_bitow(self,str):
        tab = []
        for i in str:
            tab.append(int(i))
        return tab


    probki_na_bit = int(len(dane.czas_t) / dane.ilosc_bitow)
    def kluczowanie_ASK(self,A1,A2,b):
        Za = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i/self.probki_na_bit)] == 0:
                Za.append(A1 * np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
            else:
                Za.append(A2 * np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
        print("Modulacja ASK")
        return Za

    def kluczowanie_PSK(self,b):
        Zp = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i / self.probki_na_bit)] == 0:
                Zp.append(np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i]))
            else:
                Zp.append(np.sin(2 * np.pi * self.dane.czestotliowsc_fali_nosnej_FN * self.dane.czas_t[i] + np.pi))
        return Zp

    def kluczowanie_FSK(self,b):
        Zf = []
        for i in range(len(self.dane.czas_t)):
            if b[int(i / self.probki_na_bit)] == 0:
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

    def zdemoduluj_ASK(self,zmodulowany_sygnal_bitowy_ASK):
        xt = demodulacja.pomnoz(zmodulowany_sygnal_bitowy_ASK, "ASK")
        pt = demodulacja.podziel_sygnal_i_zsumuj(xt)
        ret = demodulacja.zwroc_bitowo(pt, "ASK")
        print("Demodulacja ASK")
        return ret

    def zdemoduluj_FSK(self,sygnal_FSK):
        xt_FN1 = demodulacja.pomnoz(sygnal_FSK, "FSK")
        xt_FN2 = demodulacja.pomnoz(sygnal_FSK, "FSK", False)

        pt_FN1 = demodulacja.podziel_sygnal_i_zsumuj(xt_FN1)
        pt_FN2 = demodulacja.podziel_sygnal_i_zsumuj(xt_FN2)

        pt = demodulacja.dzialanie_na_p1_p2_FSK(pt_FN1, pt_FN2)
        ret = demodulacja.zwroc_bitowo(pt, "FSK")
        print("Demodulacja FSK")
        return ret

    def zdemoduluj_PSK(selfself,zmodulowany_sygnal_bitowy_PSK):
        xt = demodulacja.pomnoz(zmodulowany_sygnal_bitowy_PSK, "PSK")
        pt = demodulacja.podziel_sygnal_i_zsumuj(xt)
        ret = demodulacja.zwroc_bitowo(pt, "PSK")
        print("Demodulacja PSK")
        return ret

class Hamming:
    indeksy_do_usuniecia = [0, 1, 3, 7]
    c_global = []
    p = []
    def zakoduj_7_4(self, bits):
        result = np.arange(0, 7)
        if len(bits) >= 4:
            if len(bits) % 4 != 0:
                print("Bledny inpiut : Podaj wartosc podzielna przez 4")
                return -1

            kawalki = int(len(bits) / 4)
            tabs = np.array_split(bits, kawalki)
            result = []
            for i in tabs:
                z = [
                    i[0] ^ i[1] ^ i[3],
                    i[0] ^ i[2] ^ i[3],
                    i[0],
                    i[1] ^ i[2] ^ i[3],
                    i[1],
                    i[2],
                    i[3]
                ]
                result += z
        return result
    def dekoduj_7_4(self, bits):
        # print(f"Negacja:             {self.zmien_bit_dla_7_4(2,bits)}")
        result = []
        if len(bits) >= 7:
            kawalki = int(len(bits) / 7)
            tabs = np.array_split(bits, kawalki)
            for i in tabs:
                prim1 = i[2] ^ i[4] ^ i[6]
                prim2 = i[2] ^ i[5] ^ i[6]
                prim3 = i[4] ^ i[5] ^ i[6]

                dasz1 = prim1 ^ i[0]
                dasz2 = prim2 ^ i[1]
                dasz3 = prim3 ^ i[3]

                s = dasz1 * 1 + dasz2 * 2 + dasz3 * 4
                if s != 0:
                    # print("HALO")
                    i[s - 1] = int(not i[s - 1])

                result.append(i[2])
                result.append(i[4])
                result.append(i[5])
                result.append(i[6])

        return result
    def zmien_bit_dla_7_4(self,index,bits):
        if index == 1:
            bits[2] = int(not bits[2])
        elif index == 2:
            bits[4] = int(not bits[4])
        elif index == 3:
            bits[5] = int(not bits[5])
        elif index == 4:
            bits[6] = int(not bits[6])

        return bits
    def wyznacz_p(self):
        p = []
        for liczba in range(1, 16):
            liczba_bitowa = bin(liczba)[2:]
            liczba_bitowa = liczba_bitowa.zfill(4)  # Wypełnienie zerami do 4 bitów
            z = []
            for i in liczba_bitowa:
                z.insert(0, int(i))
            p.append(z)
        self.p = np.delete(p, self.indeksy_do_usuniecia, axis=0)
    def zakoduj_15_11(self,bits):
        g = self.wyznacz_macierz_generujaca_15_11()
        c = np.matmul(bits, g) % 2 #mnożenie macierzy
        self.c_global = c
        print(f"Przed zmiana:          {self.c_global}")
        k = 4
        self.c_global[k] = int(not self.c_global[k])
        print(f"Zmieniony bit w ciagu: {self.c_global} pod indexem: {k + 1}")
        return c
    def wyznacz_macierz_generujaca_15_11(self):
        self.wyznacz_p()
        i = np.eye(15 - len(self.indeksy_do_usuniecia))
        g = np.concatenate((self.p, i), axis=1)
        return g
    def wyznacz_wektor_syndromu(self):
        h = np.concatenate((np.eye(4), np.transpose(self.p)), axis=1)
        wektor_s = np.matmul(self.c_global, np.transpose(h)) % 2

        return wektor_s
    def odkoduj_15_11(self,bits):
        wektor_s = self.wyznacz_wektor_syndromu()
        s = int(wektor_s[0] * 1 + wektor_s[1] * 2 + wektor_s[2] * 4 + wektor_s[3] * 8)
        if s > 0:
            print(f"Bit {self.znajdz_index(s) - 3} jest bledny - zmiana wartosci...")
            s = self.znajdz_index(s)
            bits[s] = int(not bits[s])
            bits = np.delete(bits, [0, 1, 2, 3], axis=0)

        return bits
    def znajdz_index(self,s):
        indeksy = {
            1: 0,
            2: 1,
            3: 4,
            4: 2,
            5: 5,
            6: 6,
            7: 7,
            8: 3,
            9: 8,
            10: 9,
            11: 10,
            12: 11,
            13: 12,
            14: 13,
            15: 14
        }
        return indeksy[s]

class PomocDoLab7:
    a = [0,1,2]
    def generuj_bialy_szum(self,):
        czas = demodulacja.dane.czas_calkowity_TC
        czestotliwosc_probkowania_FS = demodulacja.dane.czestotliwosc_probkowania_FS
        ilosc_probek = int(czas * czestotliwosc_probkowania_FS)
        min = -1
        max = 1
        szum = np.random.uniform(min, max, ilosc_probek)
        return szum

    def pomnoz_gt_przez_alfa_beta(self, gt , alfa_beta):
        return np.multiply(gt, alfa_beta)

    def dodaj_szum_do_sygnalu(self, sygnal, szum, alfa):
        szum_z_alfa = self.pomnoz_gt_przez_alfa_beta(szum, alfa)
        print("Dodanie bialego szumu do sygnalu")
        return sygnal + szum_z_alfa


    def znormalizuj_sygnal(self, sygnal):
        max_val = np.max(sygnal)
        min_val = np.min(sygnal)
        znormalizowany_sygnal = 2 * (sygnal - min_val) / (max_val - min_val) - 1
        return znormalizowany_sygnal

    def licz_rozne_elementy(self, tablica1, tablica2):
        roznice = 0
        for i in range(len(tablica1)):
            if tablica1[i] != tablica2[i]:
                roznice += 1

        return roznice

    def demodulacja(self, sygnal, XSK):
        if XSK == "ASK":
            bity_po_demodulacji = demodulacja.zdemoduluj_ASK(sygnal)
        elif XSK == "FSK":
            bity_po_demodulacji = demodulacja.zdemoduluj_FSK(sygnal)
        else:
            bity_po_demodulacji = demodulacja.zdemoduluj_PSK(sygnal)

        return bity_po_demodulacji

    def modulacja_z_szumem(self,zakodowany_sygnal_bitowy, XSK, alfa):
        if XSK == "ASK":
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_ASK(0.5, 1, zakodowany_sygnal_bitowy)
        elif XSK == "FSK":
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_FSK(zakodowany_sygnal_bitowy)
        else:
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_PSK(zakodowany_sygnal_bitowy)

        bialy_szum = pomoc.generuj_bialy_szum()
        zmodulowany_sygnal_bitowy_XSK = pomoc.znormalizuj_sygnal(zmodulowany_sygnal_bitowy_XSK)
        sygnal_z_szumem = pomoc.dodaj_szum_do_sygnalu(zmodulowany_sygnal_bitowy_XSK, bialy_szum, alfa)

        return sygnal_z_szumem

    def modulacja_z_tlumieniem(self,zakodowany_sygnal_bitowy, XSK, beta):
        if XSK == "ASK":
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_ASK(0.5, 1, zakodowany_sygnal_bitowy)
        elif XSK == "FSK":
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_FSK(zakodowany_sygnal_bitowy)
        else:
            zmodulowany_sygnal_bitowy_XSK = kluczowanie.kluczowanie_PSK(zakodowany_sygnal_bitowy)

        gt = pomoc.generuj_gt(beta)
        # zmodulowany_sygnal_bitowy_XSK = pomoc.znormalizuj_sygnal(zmodulowany_sygnal_bitowy_XSK)
        sygnal_z_ukladem = pomoc.pomnoz_gt_do_sygnalu(zmodulowany_sygnal_bitowy_XSK, gt)
        return sygnal_z_ukladem

    def oblicz_wspolczynnik_BER(self,ciag_bitow_w_tab, zdekodowane):
        return round(self.licz_rozne_elementy(ciag_bitow_w_tab, zdekodowane) / len(ciag_bitow_w_tab), 2) * 100
    def generuj_gt(self, beta):
        return np.e **(-beta * demodulacja.dane.czas_t)

    def dodaj_gt_do_sygnalu(self, sygnal, szum, beta):
        szum_z_alfa = self.pomnoz_gt_przez_alfa_beta(szum, beta)
        print("Dodanie bialego szumu do sygnalu")
        return sygnal + szum_z_alfa

    def pomnoz_gt_do_sygnalu(self, sygnal, gt):
        print("Pomnozenie tlumienia ")
        return sygnal * gt

    def modyfikacja_ukladu(self, ciag_bitow_w_tab, XSK, rodzaj, alfaX = 0, betaX = 0):
        alfa = alfaX
        beta = betaX
        print(f"Ciag bitow:        {ciag_bitow_w_tab}")
        print(f"Dlugosc slowa bitowego: {len(ciag_bitow_w_tab)}")
        zakodowany_sygnal_bitowy = hamming.zakoduj_7_4(ciag_bitow_w_tab)
        print(f"zakodowane bity:   {zakodowany_sygnal_bitowy}")
        # print(f"ILOSC BITOW: [{len(zakodowany_sygnal_bitowy)}]")
        if rodzaj == "szum":
            sygnal_z_X = self.modulacja_z_szumem(zakodowany_sygnal_bitowy, XSK, alfa)
        elif rodzaj == "tlumienie":
            sygnal_z_X = self.modulacja_z_tlumieniem(zakodowany_sygnal_bitowy, XSK, beta)
            plt.plot(sygnal_z_X)
            plt.show()
        elif rodzaj == "kombinacja_I":
            sygnal_z_I = self.modulacja_z_szumem(zakodowany_sygnal_bitowy, XSK, alfa)
            gt = pomoc.generuj_gt(beta)
            sygnal_z_X = pomoc.pomnoz_gt_do_sygnalu(sygnal_z_I, gt)
        elif rodzaj == "kombinacja_II":
            sygnal_z_I = self.modulacja_z_tlumieniem(zakodowany_sygnal_bitowy, XSK, beta)
            bialy_szum = pomoc.generuj_bialy_szum()
            sygnal_z_X = pomoc.dodaj_szum_do_sygnalu(sygnal_z_I, bialy_szum, alfa)
        else:
            print("Podaj prawidlowy rodzaj : szum/tlumienie/kombinacja_I/kombinacja_II")

        bity_po_demodulacji = self.demodulacja(sygnal_z_X, XSK)
        zdekodowane = hamming.dekoduj_7_4(bity_po_demodulacji)
        print(f"zdekodowane bity:  {zdekodowane}")
        print(f"Wspolczynnik afa: {alfa}")
        print(f"Wspolczynnik beta: {beta}")
        ber = pomoc.oblicz_wspolczynnik_BER(ciag_bitow_w_tab, zdekodowane)
        print(f"Wspolczynnik BER: {ber}%")

    # def modyfikacja_kombo(self, zakodowany_sygnal_bitowy, XSK, alfaX, betaX, kombinacja):







kluczowanie = Kluczowanie()
demodulacja = Demodulacja()
hamming = Hamming()
pomoc = PomocDoLab7()

ciag_bitow_w_str = kluczowanie.ASCI_na_bity("Ostatnie zadanie")
ciag_bitow_w_tab = kluczowanie.zmien_strina_na_tablice_bitow(ciag_bitow_w_str)

print("---------------------------------------ZADANIE 1-----------------------------------------")
print(f"Ciag bitow:        {ciag_bitow_w_tab}")
print(f"Dlugosc slowa bitowego: {len(ciag_bitow_w_tab)}")
zakodowany_sygnal_bitowy = hamming.zakoduj_7_4(ciag_bitow_w_tab)
print(f"zakodowane bity:   {zakodowany_sygnal_bitowy}")
# print(f"ILOSC BITOW: [{len(zakodowany_sygnal_bitowy)}]")
zmodulowany_sygnal_bitowy_ASK = kluczowanie.kluczowanie_ASK(0.5, 1, zakodowany_sygnal_bitowy)

bity_po_demodulacji = demodulacja.zdemoduluj_ASK(zmodulowany_sygnal_bitowy_ASK)
zdekodowane = hamming.dekoduj_7_4(bity_po_demodulacji)
print(f"zdekodowane bity:  {zdekodowane}")
print(f"Ilosc rozniacych sie elementow przed i po: {pomoc.licz_rozne_elementy(ciag_bitow_w_tab, zdekodowane)}")

# alfa = 2
# print("\n\n---------------------------------------ZADANIE 2 ASK -----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "ASK", "szum", alfa)
#
# print("\n\n---------------------------------------ZADANIE 2 PSK -----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "PSK", "szum", alfa)
#
# print("\n\n---------------------------------------ZADANIE 2 FSK -----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "FSK", "szum", alfa)



beta = 20
print("\n\n---------------------------------------ZADANIE 3 ASK-----------------------------------------")

pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "ASK", "tlumienie", 0, beta)

print("\n\n---------------------------------------ZADANIE 3 PSK-----------------------------------------")

# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "PSK", "tlumienie", 0, beta)
#
# print("\n\n---------------------------------------ZADANIE 3 FSK-----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "FSK", "tlumienie", 0, beta)
#
#
# beta = 20
# alfa = 2
# print("\n\n---------------------------------------ZADANIE 4 ASK kombinacja_I-----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "ASK", "kombinacja_I", alfa, beta)
#
# print("\n\n---------------------------------------ZADANIE 4 ASK kombinacja_II-----------------------------------------")
#
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "ASK", "kombinacja_II", alfa, beta)
#
#
# print("\n\n---------------------------------------ZADANIE 4 PSK kombinacja_I-----------------------------------------")
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "PSK", "kombinacja_I", alfa, beta)
#
# print("\n\n---------------------------------------ZADANIE 4 PSK kombinacja_II-----------------------------------------")
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "PSK", "kombinacja_I", alfa, beta)
#
# print("\n\n---------------------------------------ZADANIE 4 FSK kombinacja_I-----------------------------------------")
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "FSK", "kombinacja_I", alfa, beta)
#
# print("\n\n---------------------------------------ZADANIE 4 FSK kombinacja_II-----------------------------------------")
# pomoc.modyfikacja_ukladu(ciag_bitow_w_tab, "FSK", "kombinacja_I", alfa, beta)
#


