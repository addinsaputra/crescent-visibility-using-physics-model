Bintang adalah *point source*, sehingga pembesaran teleskop meredupkan langit () tetapi tidak meredupkan bintangnya (sehingga kontras meningkat).

Untuk hilal (*extended source*), pembesaran teleskop akan **meredupkan keduanya** (baik kecerahan permukaan Bulan maupun langit latar belakang). 

Dalam konteks hilal, teleskop membantu bukan dengan membuat hilal "lebih terang" secara intrinsik per satuan luas, tetapi dengan **memperbesar ukuran sudut** hilal agar melampaui ambang batas resolusi mata dan meningkatkan pemisahan visual dari silau (*glare*). 

Hilal sangat tipis. Meskipun ia *extended* secara panjang, lebarnya sering kali di bawah resolusi mata.

* Jika lebar sabit hilal  (menit busur (arcmin2)) setelah diperbesar teleskop, mata mungkin masih menganggap lebarnya sebagai *point source* (mengikuti Hukum Ricco).
* Oleh karena itu, perhitungan visibilitas hilal harus memperhitungkan **lebar sabit (W)**. Jika  masih sangat kecil, rumus *point source* sebagian masih berlaku untuk dimensi lebarnya.

import math

class TelescopeVisibilityModel:
    def __init__(self):
        """
        Implementasi faktor koreksi visibilitas teleskop berdasarkan 
        Paper: TELESCOPIC LIMITING MAGNITUDES (Schaefer, 1990).
        """
        # [cite_start]1. Fb: Faktor Binokular (Binocular Factor) [cite: 267]
        # Sensitivitas 2 mata vs 1 mata (teleskop monokular).
        self.Fb = 1.41

    def calculate_factors(self, D, Ds, M, De, t1, n, theta):
        """
        Menghitung variabel-variabel koreksi optik.
        
        Parameters:
        -----------
        D     : Aperture/Diameter teleskop 
        Ds    : Diameter halangan cermin sekunder
        M     : Magnifikasi / Perbesaran
        De    : Diameter pupil mata pengamat = 7 exp (-0.5(A/100)^2)
        A     : usia pengamat (dalam tahun)
        t1    : Transmisi per permukaan optik
        n     : Jumlah permukaan optik (lensa/cermin)
        theta : Ukuran seeing disk (dalam detik busur / arcseconds)
        
        Returns:
        --------
        Dictionary berisi nilai-nilai variabel yang diminta.
        """
        
        # [cite_start]2. (Ds/D)^2: Fraksi Kehilangan Cahaya Obstruksi [cite: 282]
        # Cahaya yang terhalang oleh cermin sekunder.
        obstruction_fraction = (Ds / D) ** 2
        
        # [cite_start]3. Ft: Faktor Transmisi Teleskop [cite: 287]
        # Memperhitungkan kehilangan cahaya pada setiap permukaan (t1^n) 
        # dan halangan sekunder (1 - obstruction_fraction).
        # Rumus: Ft = 1 / (t1^n * (1 - (Ds/D)^2)) #Fraksi Kehilangan Cahaya Obstruksi (hanya berlaku untuk teleskop reflektor yang mengandalkan cermin sekunder, untuk teleskop refraktor, Ds/D = 0)
        try:
            Ft = 1 / ((t1 ** n) * (1 - obstruction_fraction))
        except ZeroDivisionError:
            Ft = float('inf') # Menghindari error jika transmisi 0
            
        # [cite_start]4. Fp: Faktor Pupil [cite: 300-302]
        # Koreksi jika exit pupil teleskop (D/M) lebih besar dari pupil mata (De).
        # Rumus: Fp = (D / (M * De))^2 jika De < D/M, else 1.0
        exit_pupil_telescope = D / M
        if De < exit_pupil_telescope:
            Fp = (D / (M * De)) ** 2
        else:
            Fp = 1.0
            
        # [cite_start]5. Fa: Faktor Bukaan (Aperture Factor) [cite: 307]
        # Rasio pengumpulan cahaya mata vs teleskop.
        # Rumus: Fa = (De / D)^2
        Fa = (De / D) ** 2
        
        # [cite_start]6. Fm: Faktor Magnifikasi Latar Belakang [cite: 309]
        # Pengurangan brightness langit akibat perbesaran.
        # Rumus: Fm = 1/M^2
        Fm = 1/M**2
        
        # [cite_start]7. Fr: Faktor Resolusi Mata [cite: 319-320]
        # Koreksi jika objek titik diperbesar hingga tampak buram/besar 
        # (melebihi batas resolusi mata ~900 detik busur).
        # Rumus: Fr = sqrt(2 * theta * M / 900) jika > 900", else 1.0
        apparent_size = 2 * theta * M
        if apparent_size > 900:
            Fr = (apparent_size / 900) ** 0.5
        else:
            Fr = 1.0
            
        return {
            "Fb (Binocular Factor)": self.Fb,
            "Obstruction Fraction (Ds/D)^2": obstruction_fraction,
            "t1 (Transmission per surface)": t1,
            "n (Number of surfaces)": n,
            "Ft (Telescope Transmission Factor)": Ft,
            "Fp (Pupil Factor)": Fp,
            "Fa (Aperture Factor)": Fa,
            "Fm (Magnification Factor)": Fm,
            "Fr (Resolution Factor)": Fr
        }

# --- Contoh Penggunaan ---
if __name__ == "__main__":
    # Inisialisasi model
    model = TelescopeVisibilityModel()
    
    # Contoh parameter menggunakan teleskop refraktor
    D_val = 150.0      # Diameter teleskop 150 mm
    Ds_val = 0.0       # Diameter sekunder 0 mm
    M_val = 50.0       # Perbesaran 50x
    De_val = 7.0       # Pupil mata 7 mm (adaptasi gelap)
    t1_val = 0.95      # Transmisi per lensa/cermin 95%
    n_val = 6          # 2 cermin + 4 elemen lensa eyepiece
    theta_val = 2.0    # Seeing 2 detik busur
    
    # Hitung faktor
    results = model.calculate_factors(D_val, Ds_val, M_val, De_val, t1_val, n_val, theta_val)
    
    # Tampilkan hasil
    print("Hasil Perhitungan Variabel Schaefer (1990):")
    for key, value in results.items():
        print(f"{key}: {value:.4f}")



rumus koreksi teleskop adalah:
FB = Fb * Ft * Fp * Fa * Fm #ini adalah rumus koreksi teleskop untuk sky brightness 
FI = Fb * Ft * Fp * Fa * Fr #ini adalah rumus koreksi teleskop untuk luminansi hilal

B_eff = B_0 / FB #ini adalah rumus koreksi teleskop final untuk sky brightness
I_eff = I_0 / FI #ini adalah rumus koreksi teleskop final untuk luminansi hilal

visibilitas hilal dengan koreksi teleskop adalah:
delta m = 2.5 log (I_eff/B_eff) #ini adalah rumus koreksi teleskop final untuk visibilitas hilal
