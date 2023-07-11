import numpy as np
import matplotlib.pyplot as plt
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
        print(f"Negacja:             {self.zmien_bit_dla_7_4(2,bits)}")
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

x = Hamming()
bits = [1, 0, 1, 1, 1 ,0 ,1 ,0]
print(f"Bity do zakodowania: {bits}")
print("-------")
zakodowane = x.zakoduj_7_4(bits)
print(f"zakodowane bity    : {zakodowane}")
zdekodowane = x.dekoduj_7_4(zakodowane)
print(f"zdekodowane bity   : {zdekodowane}")

print("-------")
print("-------")
print("-------")


bits = [0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1]
print(f"bity do zakodowania: {bits}")
print("-------")
zakodowne = x.zakoduj_15_11(bits)
print(f"Zakodowane bity:       {zakodowne}")
odkodowane = x.odkoduj_15_11(zakodowne)
print(f"Odkodowane bity:       {odkodowane}")

