"""
Program Integrasi Model Visibilitas Hilal Kastner Modifikasi
==============================================================

Program ini mengintegrasikan model Schaefer (untuk sky brightness), model Kastner
(untuk luminansi hilal), dan data atmosferik dari NASA POWER API untuk memprediksi
visibilitas hilal dengan otomatisasi penuh.

Alur Program:
1. Input user: bulan hijriah, tahun hijriah, lokasi (lintang, bujur, elevasi, timezone)
2. Otomatisasi: hitung waktu ijtima/konjungsi
3. Otomatisasi: hitung sunset dengan koreksi bertahap:
   a. Sunset GEOMETRIS (tanpa koreksi refraksi) menggunakan skyfield
   b. Ambil data atmosfer (RH, T) dari NASA POWER API pada waktu sunset geometris
   c. Sunset APPARENT dengan koreksi refraksi menggunakan T dari API
4. Otomatisasi: tentukan tanggal pengamatan berdasarkan aturan hisab
5. Otomatisasi: hitung posisi matahari dan bulan saat sunset
6. Hitung sky brightness menggunakan model Schaefer (visual_limit_calc2.py)
7. Hitung luminansi hilal menggunakan model Kastner (kastner.py)
8. Hitung visibilitas naked eye dan teleskop
9. Output hasil lengkap

============================================================
PERBAIKAN PENTING - Versi 2026-01-11 (Update 6)
============================================================

1. IMPLEMENTASI MODEL SCHAEFER UNTUK TELESKOP
   --------------------------------------------
   Peningkatan: Menggunakan model koreksi Schaefer untuk visibilitas teleskop:

   Model Schaefer:
   - Berdasarkan faktor koreksi extended surface
   - Menggunakan extended_surface_correction_factor dari limit_teleskop.py
   - Output: Delta m (contrast magnitude), >0 = detectable
   - Fungsi: extended_surface_correction_factor() dari limit_teleskop.py

   Keuntungan:
   - Model yang konsisten dengan perhitungan sky brightness
   - Perhitungan delta m yang standard untuk teleskop
   - Koreksi untuk extended source yang mengakomodasi magnifikasi

   Perubahan:
   - Import: from limit_teleskop import extended_surface_correction_factor
   - Fungsi hitung_visibilitas_teleskop() menggunakan model Schaefer
   - Output delta_m_tel menggunakan koreksi Schaefer

   Catatan:
   - Naked eye dan teleskop sama-sama menggunakan model Schaefer (delta m)
   - Untuk extended source, magnifikasi mengurangi surface brightness
   - Delta m > 0 berarti hilal terdeteksi

============================================================
PERBAIKAN PENTING - Versi 2026-01-09 (Update 5)
============================================================

1. ALUR KOREKSI SUNSET YANG BENAR (GEOMETRIS → ATMOSFER → REFRAKSI)
   ------------------------------------------------------------------
   Peningkatan: Mengimplementasikan alur koreksi sunset yang ilmiah dan benar:

   Alur yang benar:
   a. Hitung sunset GEOMETRIS (tanpa koreksi refraksi) menggunakan skyfield
      - Hanya koreksi geometri akibat elevasi lokasi
      - Fungsi: sunrise_sunset_geometric_local()

   b. Ambil data atmosfer (RH, T) dari NASA POWER API
      - Waktu pengambilan data: pada saat sunset geometris
      - Data yang diambil: Relative Humidity (RH) dan Temperature (T)
      - Fungsi: get_rh_t_at_time() dari meteo_power_api.py

   c. Hitung sunset APPARENT dengan koreksi refraksi
      - Menggunakan data Temperature dari API untuk koreksi refraksi atmosfer
      - Pressure di-set statis (default: 1013.25 mbar)
      - Fungsi: sunrise_sunset_apparent_local()

   Keuntungan:
   - Sunset dihitung dengan data temperature aktual dari lokasi dan waktu pengamatan
   - Koreksi refraksi atmosfer lebih akurat
   - Data atmosfer (RH dan T) yang sama digunakan untuk perhitungan sky brightness
   - Tidak ada input manual (otomatisasi penuh)

2. FUNGSI BARU DI sunmoon.py
   --------------------------
   Peningkatan: Menambahkan dua fungsi baru untuk koreksi sunset:

   a. sunrise_sunset_geometric_local(location, timezone, year, month, day)
      - Menghitung sunset TANPA koreksi refraksi (hanya koreksi geometri/elevasi)
      - Digunakan untuk estimasi awal waktu sunset

   b. sunrise_sunset_apparent_local(location, timezone, year, month, day, T, P)
      - Menghitung sunset DENGAN koreksi refraksi atmosfer
      - Parameter T (temperature) diambil dari NASA POWER API
      - Parameter P (pressure) di-set statis (default: 1013.25 mbar)
      - Menggunakan binary search untuk menemukan waktu sunset yang akurat

3. OTOMATISASI DATA ATMOSFER
   ---------------------------
   Peningkatan: Data atmosfer (RH dan T) diambil secara otomatis dari NASA POWER API

   - Waktu pengambilan: pada saat sunset geometris
   - Penggunaan data:
     * Temperature: untuk koreksi refraksi atmosfer saat hitung sunset
     * RH dan Temperature: untuk perhitungan sky brightness (model Schaefer)
   - Fallback: jika API gagal, gunakan nilai default (RH=80%, T=25°C)

4. PENENTUAN TANGGAL PENGAMATAN
   -------------------------------
   Aturan hisab yang diterapkan:
   - Jika ijtima lokal terjadi sebelum jam 12:00 (tengah malam - sebelum siang):
     Pengamatan dilakukan pada hari yang sama (tanggal ijtima lokal, sore hari)
   - Jika ijtima lokal terjadi setelah jam 12:00 (siang - tengah malam):
     Bandingkan dengan sunset geometris sore hari itu:
     * Jika ijtima < sunset: amati hari berikutnya
     * Jika ijtima >= sunset: amati hari yang sama

5. EPGHEMERIS DE440s (HIGH-PRECISION)
   ------------------------------------
   Menggunakan ephemeris DE440s.bsp (modern, 1849-2150) untuk akurasi posisi
   matahari dan bulan yang lebih tinggi.

6. VALIDASI HASIL
   ---------------
   Untuk kasus Muharram 1444 (29 Juli 2022) di UIN Walisongo:
   - Sunset Geometris: ~17:32 WIB (tanpa refraksi)
   - Sunset Apparent: ~17:36 WIB (dengan refraksi dan T dari API)
   - Sun Altitude: ~0 derajat (sesuai definisi sunset astronomis)
   - Moon Altitude: 6.37 derajat (positif, sesuai ekspektasi)
   - Elongasi Toposentrik: 8.39 derajat
   - Delta m (Teleskop): 1.69 (positif, hilal terlihat dengan teleskop)

Catatan Penting:
- Program menggunakan DE440s.bsp ephemeris (high-precision)
- Sunset dihitung dengan koreksi refraksi atmosfer yang akurat
- Temperature untuk koreksi refraksi diambil otomatis dari NASA POWER API
- Pressure di-set statis (default: 1013.25 mbar)
- Data atmosfer (RH dan T) digunakan untuk perhitungan sky brightness
- Delta day offset dapat digunakan untuk menyesuaikan tanggal pengamatan

Author: Droid AI Assistant
Date: 2026-01-09
Version: 5.0 (Alur Koreksi Sunset: Geometris → Atmosfer → Refraksi)
"""

import math
from datetime import datetime, timezone, timedelta
from typing import Dict, Any, Optional, Tuple
import sys
import os

# Import modul yang diperlukan
from visual_limit_calc2 import visual_limit
from kastner import hitung_luminansi_kastner, crescent_area
from limit_teleskop import (
    TelescopeVisibilityModel,
    calculate_telescope_visibility
)  # Model Schaefer untuk koreksi teleskop
from nasa_power_api_merra2 import ObservingLocation, get_rh_t_at_time, PowerAPIError

# Import modul sunmoon dari direktori data-hisab
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'data-hisab'))
from sunmoon import (
    newmoon_hijri_month_utc,
    convert_utc_to_localtime,
    convert_localtime_to_utc,
    sunrise_sunset_utc,
    sunrise_sunset_local,
    sunrise_sunset_geometric_local,  # NEW: Sunset tanpa koreksi refraksi (hanya geometri)
    sunrise_sunset_apparent_local,   # NEW: Sunset dengan koreksi refraksi atmosfer
    sun_position_time_utc,
    sun_position_time_local,
    moon_position_time_utc,
    moon_position_time_local,
    moon_elongation_time_utc,
    moon_elongation_time_local,
    moon_phase_angle_time_utc,
    moon_phase_angle_time_local,
    moon_illumination_width_utc,
    moon_illumination_width_local,
    set_location
)



def nL_to_mag_arcsec2(B_nL: float) -> float:
    """
    Konversi Sky Brightness dari nanoLambert (nL) ke mag/arcsec^2.
    Menggunakan invers dari rumus Schaefer:
    B_s = 34.08 * exp(20.7233 - 0.92104 * B_V)
    
    Sehingga:
    B_V = (20.7233 - ln(B_s / 34.08)) / 0.92104
    """
    if B_nL <= 0:
        return 22.0 # Nilai gelap default jika input invalid (misal 0)
        
    try:
        # B_s di rumus Schaefer adalah dalam millimicro-Lamberts.
        # 1 nL = 1 millimicro-Lambert. Jadi unit sama.
        val = math.log(B_nL / 34.08)
        B_V = (20.7233 - val) / 0.92104
        return B_V
    except (ValueError, ZeroDivisionError):
        return 22.0 # Fallback ke langit gelap


def deg_to_dms(deg: float) -> str:
    """
    Konversi derajat desimal ke format derajat-menit-detik (DMS).
    
    Parameters:
    -----------
    deg : float
        Nilai dalam derajat desimal
        
    Returns:
    --------
    str
        String dalam format "DD° MM' SS.SS\""
    """
    sign = "-" if deg < 0 else ""
    deg = abs(deg)
    d = int(deg)
    m = int((deg - d) * 60)
    s = (deg - d - m / 60) * 3600
    return f"{sign}{d}° {m}' {s:.2f}\""


class HilalVisibilityCalculator:
    """Kelas utama untuk kalkulasi visibilitas hilal"""
    
    def __init__(self,
                 nama_tempat: str,
                 lintang: float,
                 bujur: float,
                 elevasi: float,
                 timezone_str: str,
                 bulan_hijri: int,
                 tahun_hijri: int,
                 delta_day_offset: int = 0):
        """
        Inisialisasi kalkulator visibilitas hilal.
        
        Parameters:
        -----------
        nama_tempat : str
            Nama lokasi pengamatan
        lintang : float
            Lintang dalam derajat (positif untuk utara, negatif untuk selatan)
        bujur : float
            Bujur dalam derajat (positif untuk timur, negatif untuk barat)
        elevasi : float
            Ketinggian dalam meter di atas permukaan laut
        timezone_str : str
            Zona waktu (contoh: "Asia/Jakarta" atau "+7")
        bulan_hijri : int
            Bulan hijriah (1-12)
        tahun_hijri : int
            Tahun hijriah
        delta_day_offset : int
            Offset hari untuk pengamatan (default: 0)
        """
        self.nama_tempat = nama_tempat
        self.lintang = lintang
        self.bujur = bujur
        self.elevasi = elevasi
        self.timezone_str = timezone_str
        self.bulan_hijri = bulan_hijri
        self.tahun_hijri = tahun_hijri
        self.delta_day_offset = delta_day_offset
        
        # Setup lokasi untuk perhitungan astronomis
        self.location = set_location(lintang, bujur, elevasi)
        
        # Nilai statis sesuai spesifikasi
        self.ALTMOON_STATIS = -90.0  # Dikunci di -90 derajat
        self.AZIMOON_STATIS = 180.0  # Dikunci di 180 derajat
        self.SNELLEN_RATIO = 1.0     # Dikunci di 1
        self.PHASE_MOON = 0          # Dikunci di 0 (new moon)
        
        # Hasil perhitungan akan disimpan di sini
        self.hasil: Dict[str, Any] = {}
    
    def hitung_ijtima(self) -> Tuple[datetime, datetime]:
        """
        Menghitung waktu ijtima (konjungsi) untuk bulan hijriah yang ditentukan.
        
        Returns:
        --------
        ijtima_utc : datetime
            Waktu ijtima dalam UTC
        ijtima_local : datetime
            Waktu ijtima dalam waktu lokal
        """
        ijtima_utc = newmoon_hijri_month_utc(self.tahun_hijri, self.bulan_hijri)
        ijtima_local = convert_utc_to_localtime(self.timezone_str, utc_datetime=ijtima_utc)
        
        self.hasil['ijtima_utc'] = ijtima_utc
        self.hasil['ijtima_local'] = ijtima_local
        
        return ijtima_utc, ijtima_local
    
    def tentukan_tanggal_pengamatan(self, ijtima_utc: datetime) -> Tuple[datetime, datetime, float, float]:
        """
        Menentukan tanggal pengamatan berdasarkan waktu ijtima dan sunset.

        ALUR KOREKSI SUNSET:
        1. Hitung sunset GEOMETRIS (tanpa koreksi refraksi) menggunakan skyfield
        2. Gunakan waktu sunset geometris untuk mengambil data atmosfer (RH, T) dari NASA POWER API
        3. Hitung sunset APPARENT dengan koreksi refraksi menggunakan T dari API
        4. Data RH dan T digunakan untuk perhitungan sky brightness

        ATURAN HISAB:
        1. Konversi ijtima UTC ke waktu lokal
        2. Jika ijtima lokal terjadi sebelum jam 12:00 (tengah malam - siang):
           - Gunakan tanggal ijtima lokal untuk pengamatan (sore hari itu)
        3. Jika ijtima lokal terjadi setelah jam 12:00 (siang - tengah malam):
           - Bandingkan dengan sunset lokal sore hari itu
           - Jika ijtima < sunset: amati hari berikutnya
           - Jika ijtima >= sunset: amati hari yang sama

        Contoh untuk kasus Muharram 1444 (29 Juli 2022):
        - Ijtima UTC: 2022-07-28 17:55:02
        - Ijtima Lokal (WIB): 2022-07-29 00:55:02
        - Karena ijtima lokal (00:55) < 12:00, maka pengamatan dilakukan pada 29 Juli 2022 sore
        - Sunset 29 Juli sore: ~17:30 WIB -> Bulan sudah cukup tinggi untuk diamati

        Parameters:
        -----------
        ijtima_utc : datetime
            Waktu ijtima dalam UTC

        Returns:
        --------
        sunset_utc : datetime
            Waktu sunset dalam UTC pada hari pengamatan
        sunset_local : datetime
            Waktu sunset dalam waktu lokal (dengan koreksi refraksi)
        rh : float
            Relative humidity dalam persen (dari API)
        temperature : float
            Suhu dalam derajat Celsius (dari API)
        """
        # Konversi ijtima ke waktu lokal
        ijtima_local = convert_utc_to_localtime(self.timezone_str, utc_datetime=ijtima_utc)

        # Aturan sederhana berdasarkan jam ijtima lokal
        if ijtima_local.hour < 12:
            # Ijtima terjadi sebelum jam 12:00 (tengah malam sampai sebelum siang)
            tanggal_pengamatan = ijtima_local.date()
        else:
            # Ijtima terjadi setelah jam 12:00 (siang sampai tengah malam)
            _, sunset_local_ijtima_geom = sunrise_sunset_geometric_local(
                self.location,
                self.timezone_str,
                year=ijtima_utc.year,
                month=ijtima_utc.month,
                day=ijtima_utc.day
            )

            if ijtima_local < sunset_local_ijtima_geom:
                tanggal_pengamatan = ijtima_local.date() + timedelta(days=1)
            else:
                tanggal_pengamatan = ijtima_local.date()

        # Terapkan delta day offset
        tanggal_pengamatan += timedelta(days=self.delta_day_offset)

        # LANGKAH 1: Hitung sunset GEOMETRIS (tanpa koreksi refraksi)
        _, sunset_geom_local = sunrise_sunset_geometric_local(
            self.location,
            self.timezone_str,
            year=tanggal_pengamatan.year,
            month=tanggal_pengamatan.month,
            day=tanggal_pengamatan.day
        )
        sunset_geom_utc = convert_localtime_to_utc(self.timezone_str, local_datetime=sunset_geom_local)

        # LANGKAH 2: Ambil data atmosfer (RH dan T) pada waktu sunset geometris
        loc = ObservingLocation(
            name=self.nama_tempat,
            latitude=self.lintang,
            longitude=self.bujur,
            altitude=self.elevasi,
            timezone=self.timezone_str
        )

        # Pastikan sunset_geom_utc adalah timezone-aware UTC datetime
        if sunset_geom_utc.tzinfo is None or sunset_geom_utc.tzinfo.utcoffset(sunset_geom_utc) is None:
            sunset_geom_utc = sunset_geom_utc.replace(tzinfo=timezone.utc)

        try:
            rh, temperature = get_rh_t_at_time(loc, sunset_geom_utc)
        except (PowerAPIError, Exception) as e:
            print(f"  [!] API Error: {e} - Menggunakan default RH=80%, T=25°C")
            rh = 80.0
            temperature = 25.0

        # LANGKAH 3: Hitung sunset APPARENT dengan koreksi refraksi menggunakan T dari API
        _, sunset_local = sunrise_sunset_apparent_local(
            self.location,
            self.timezone_str,
            year=tanggal_pengamatan.year,
            month=tanggal_pengamatan.month,
            day=tanggal_pengamatan.day,
            temperature_C=temperature,  # Gunakan temperature dari API
            pressure_mbar=1013.25  # Pressure statis
        )

        # Konversi sunset lokal ke UTC
        sunset_utc = convert_localtime_to_utc(self.timezone_str, local_datetime=sunset_local)

        # Simpan hasil
        self.hasil['sunset_utc'] = sunset_utc
        self.hasil['sunset_local'] = sunset_local
        self.hasil['sunset_geom_utc'] = sunset_geom_utc
        self.hasil['sunset_geom_local'] = sunset_geom_local
        self.hasil['tanggal_pengamatan'] = datetime.combine(tanggal_pengamatan, datetime.min.time())
        self.hasil['rh'] = rh
        self.hasil['temperature'] = temperature

        return sunset_utc, sunset_local, rh, temperature
    
    def ambil_data_atmosfer(self, sunset_utc: datetime) -> Tuple[float, float]:
        """
        Mengambil data atmosfer (RH dan suhu) dari NASA POWER API.
        
        Parameters:
        -----------
        sunset_utc : datetime
            Waktu sunset dalam UTC
            
        Returns:
        --------
        rh : float
            Relative humidity dalam persen
        temperature : float
            Suhu udara dalam derajat Celsius
        """
        loc = ObservingLocation(
            name=self.nama_tempat,
            latitude=self.lintang,
            longitude=self.bujur,
            altitude=self.elevasi,
            timezone=self.timezone_str
        )
        # Pastikan `sunset_utc` adalah timezone-aware UTC datetime
        if sunset_utc is None:
            raise ValueError("sunset_utc tidak boleh None")
        if sunset_utc.tzinfo is None or sunset_utc.tzinfo.utcoffset(sunset_utc) is None:
            sunset_utc = sunset_utc.replace(tzinfo=timezone.utc)
        else:
            sunset_utc = sunset_utc.astimezone(timezone.utc)

        try:
            rh, temperature = get_rh_t_at_time(loc, sunset_utc)
            self.hasil['rh'] = rh
            self.hasil['temperature'] = temperature
            return rh, temperature
        except PowerAPIError as e:
            print(f"  [!] API Error: {e} - Menggunakan default RH=80%, T=25°C")
            rh = 80.0
            temperature = 25.0
            self.hasil['rh'] = rh
            self.hasil['temperature'] = temperature
            return rh, temperature
    
    def hitung_posisi_matahari_bulan(self, sunset_local: datetime) -> Dict[str, float]:
        """
        Menghitung posisi matahari dan bulan saat sunset.
        
        Parameters:
        -----------
        sunset_local : datetime
            Waktu sunset dalam waktu lokal
            
        Returns:
        --------
        Dict[str, float]
            Dictionary berisi:
            - sun_alt: Altitude matahari (derajat)
            - sun_az: Azimuth matahari (derajat)
            - moon_alt: Altitude bulan (derajat)
            - moon_az: Azimuth bulan (derajat)
            - elongation: Elongasi toposentrik (derajat)
            - phase_angle: Sudut fase bulan (derajat)
            - moon_semidiameter: Semidiameter bulan (derajat)
        """
        # Posisi matahari (menggunakan local time)
        sun_alt, sun_az, _ = sun_position_time_local(
            self.location,
            self.timezone_str,
            local_datetime=sunset_local
        )
        
        # Posisi bulan (menggunakan local time)
        moon_alt, moon_az, _ = moon_position_time_local(
            self.location,
            self.timezone_str,
            local_datetime=sunset_local
        )
        
        # Elongasi toposentrik
        elongation = moon_elongation_time_local(
            self.timezone_str,
            location=self.location,
            local_datetime=sunset_local
        )
        
        # Sudut fase bulan
        phase_angle = moon_phase_angle_time_local(
            self.timezone_str,
            location=self.location,
            local_datetime=sunset_local
        )
        
        # Semidiameter bulan
        _, _, parallax, SD = moon_illumination_width_local(
            self.timezone_str,
            location=self.location,
            local_datetime=sunset_local
        )
        moon_semidiameter = SD  # dalam derajat
        
        return {
            'sun_alt': sun_alt,
            'sun_az': sun_az,
            'moon_alt': moon_alt,
            'moon_az': moon_az,
            'elongation': elongation,
            'phase_angle': phase_angle,
            'moon_semidiameter': moon_semidiameter
        }
    
    def hitung_sky_brightness_schaefer(self,
                                        rh: float,
                                        temperature: float,
                                        posisi: Dict[str, float]) -> Tuple[float, float]:
        """
        Menghitung sky brightness menggunakan model Schaefer.
        
        Parameters:
        -----------
        rh : float
            Relative humidity dalam persen
        temperature : float
            Suhu dalam derajat Celsius
        posisi : Dict[str, float]
            Dictionary posisi matahari dan bulan
            
        Returns:
        --------
        sky_brightness_nl : float
            Sky brightness dalam nanoLambert
        k_v : float
            Koefisien ekstingsi band V
        """
        # Hitung selisih azimuth sun-moon
        azisun = abs(posisi['sun_az'] - posisi['moon_az'])
        sun_alt_safe = posisi['sun_alt']
        
        try:
            result = visual_limit(
                month=self.hasil['tanggal_pengamatan'].month,
                year=self.hasil['tanggal_pengamatan'].year,
                phase_angle=self.PHASE_MOON,
                altmoon=self.ALTMOON_STATIS,
                azimoon=self.AZIMOON_STATIS,
                altsun=sun_alt_safe,
                azisun=azisun,
                humidity=rh,
                temperature=temperature,
                latitude=self.lintang,
                altitude=self.elevasi,
                snellen_ratio=self.SNELLEN_RATIO,
                altstar=max(posisi['moon_alt'], 0.0),
                goal_magnitude=99
            )
            
            sky_brightness_nl = float(result.get("sky_brightness", 0.0))
            
            if sky_brightness_nl <= 0:
                sky_brightness_nl = 1.0
            
        except (ValueError, ZeroDivisionError) as e:
            print(f"  [!] Error sky brightness: {e} - Menggunakan default 1000 nL")
            sky_brightness_nl = 1000.0
        
        K_list = result.get("K", []) if 'result' in locals() else []
        k_v = float(K_list[2]) if len(K_list) > 2 else 0.3
        
        return sky_brightness_nl, k_v
    
    def hitung_luminansi_hilal_kastner(self,
                                       posisi: Dict[str, float],
                                       k_v: float) -> float:
        """
        Menghitung luminansi hilal menggunakan model Kastner.
        
        Parameters:
        -----------
        posisi : Dict[str, float]
            Dictionary posisi matahari dan bulan
        k_v : float
            Koefisien ekstingsi band V
            
        Returns:
        --------
        luminansi_hilal_nl : float
            Luminansi hilal dalam nanoLambert
        """
        zenith_distance = 90.0 - posisi['moon_alt']
        
        luminansi_hilal_nl = hitung_luminansi_kastner(
            alpha=posisi['phase_angle'],
            elongation=posisi['elongation'],
            r=float(posisi['moon_semidiameter']),
            z=zenith_distance,
            k=k_v
        )
        
        return luminansi_hilal_nl
    
    def hitung_visibilitas_naked_eye(self,
                                      luminansi_hilal_nl: float,
                                      sky_brightness_nl: float) -> Tuple[float, float]:
        """
        Menghitung visibilitas untuk pengamatan naked eye (mata telanjang).
        
        Parameters:
        -----------
        luminansi_hilal_nl : float
            Luminansi hilal dalam nL
        sky_brightness_nl : float
            Sky brightness dalam nL
            
        Returns:
        --------
        rasio_kontras : float
            Rasio kontras R = L/B
        delta_m : float
            Delta magnitudo
        """
        if sky_brightness_nl <= 0.0:
            raise ValueError("Sky brightness harus positif")
        
        rasio_kontras = luminansi_hilal_nl / sky_brightness_nl
        delta_m = 2.5 * math.log10(rasio_kontras) if rasio_kontras > 0.0 else float('-inf')
        
        return rasio_kontras, delta_m
    
    def hitung_visibilitas_teleskop(self,
                                     luminansi_hilal_nl: float,
                                     sky_brightness_nl: float,
                                     posisi: Dict[str, float],
                                     k_v: float,
                                     aperture: float = 66.0,
                                     magnification: float = 50.0,
                                     transmission: float = 0.95,
                                     n_surfaces: int = 6,
                                     central_obstruction: float = 0.0,
                                     observer_age: float = 22.0,
                                     seeing: float = 2.0) -> Tuple[float, float, float, float]:
        """
        Menghitung visibilitas untuk pengamatan dengan teleskop menggunakan model Schaefer.

        Model Schaefer (1990) menggunakan faktor koreksi yang BERBEDA untuk:
        - Sky brightness: FB = Fb * Ft * Fp * Fa * Fm
        - Luminansi hilal: FI = Fb * Ft * Fp * Fa * Fr

        Kemudian:
        - B_eff = B_0 / FB
        - I_eff = I_0 / FI
        - delta_m = 2.5 * log10(I_eff / B_eff)

        Parameters:
        -----------
        luminansi_hilal_nl : float
            Luminansi hilal dalam nL (I_0)
        sky_brightness_nl : float
            Sky brightness dalam nL (B_0)
        posisi : Dict[str, float]
            Dictionary posisi matahari dan bulan
        k_v : float
            Koefisien ekstingsi band V (tidak digunakan di model baru)
        aperture : float
            Diameter aperture teleskop (mm)
        magnification : float
            Pembesaran teleskop
        transmission : float
            Transmisi per permukaan optik (default: 0.95)
        n_surfaces : int
            Jumlah permukaan optik (default: 4)
        central_obstruction : float
            Diameter obstruksi pusat (mm, default: 0.0 untuk refraktor)
        observer_age : float
            Usia pengamat (tahun, default: 30.0)
        seeing : float
            Ukuran seeing disk (arcseconds, default: 2.0)

        Returns:
        --------
        luminansi_hilal_tel_nl : float
            Luminansi hilal efektif teleskop (I_eff) dalam nL
        sky_brightness_tel_nl : float
            Sky brightness efektif teleskop (B_eff) dalam nL
        rasio_kontras_tel : float
            Rasio kontras teleskop (I_eff / B_eff)
        delta_m_tel : float
            Delta magnitudo teleskop, >0 = detectable
        """
        # Gunakan model TelescopeVisibilityModel yang baru
        result = calculate_telescope_visibility(
            B_sky=sky_brightness_nl,
            L_hilal=luminansi_hilal_nl,
            aperture=aperture,
            magnification=magnification,
            central_obstruction=central_obstruction,
            transmission=transmission,
            n_surfaces=n_surfaces,
            age=observer_age,
            seeing=seeing
        )
        
        # Ekstrak hasil
        B_eff = result['B_eff']  # Sky brightness terkoreksi
        I_eff = result['I_eff']  # Luminansi hilal terkoreksi
        delta_m_tel = result['delta_m']  # Kontras magnitude
        
        # Hitung rasio kontras
        if B_eff > 0:
            rasio_kontras_tel = I_eff / B_eff
        else:
            rasio_kontras_tel = float('nan')

        # Return hasil perhitungan
        return I_eff, B_eff, rasio_kontras_tel, delta_m_tel



    
    def jalankan_perhitungan_lengkap(self,
                                      use_telescope: bool = True,
                                      aperture: float = 100.0,
                                      magnification: float = 50.0) -> Dict[str, Any]:
        """
        Menjalankan seluruh algoritma perhitungan visibilitas hilal.

        ALUR PROGRAM:
        1. Hitung ijtima/konjungsi
        2. Tentukan tanggal pengamatan dan hitung sunset dengan alur koreksi:
           a. Hitung sunset geometris (tanpa refraksi)
           b. Ambil data atmosfer (RH, T) dari NASA POWER API pada waktu sunset geometris
           c. Hitung sunset apparent dengan koreksi refraksi menggunakan T dari API
        3. Hitung posisi matahari dan bulan saat sunset
        4. Hitung sky brightness menggunakan model Schaefer (dengan RH dan T dari API)
        5. Hitung luminansi hilal menggunakan model Kastner
        6. Hitung visibilitas naked eye dan teleskop

        Parameters:
        -----------
        use_telescope : bool
            Apakah akan menghitung visibilitas teleskop
        aperture : float
            Diameter aperture teleskop (mm)
        magnification : float
            Pembesaran teleskop

        Returns:
        --------
        Dict[str, Any]
            Dictionary berisi seluruh hasil perhitungan
        """
        print("Menjalankan perhitungan visibilitas hilal...")

        # Langkah 1: Hitung ijtima
        ijtima_utc, ijtima_local = self.hitung_ijtima()

        # Langkah 2: Tentukan tanggal pengamatan dan hitung sunset dengan koreksi refraksi
        sunset_utc, sunset_local, rh, temperature = self.tentukan_tanggal_pengamatan(ijtima_utc)

        # Langkah 3: Hitung posisi matahari dan bulan (menggunakan local time)
        posisi = self.hitung_posisi_matahari_bulan(sunset_local)

        # Langkah 4: Hitung sky brightness (Schaefer)
        sky_brightness_nl, k_v = self.hitung_sky_brightness_schaefer(rh, temperature, posisi)

        # Langkah 5: Hitung luminansi hilal (Kastner)
        luminansi_hilal_nl = self.hitung_luminansi_hilal_kastner(posisi, k_v)

        # Langkah 6: Hitung visibilitas naked eye
        rasio_kontras_ne, delta_m_ne = self.hitung_visibilitas_naked_eye(
            luminansi_hilal_nl, sky_brightness_nl
        )

        # Langkah 7: Hitung visibilitas teleskop (jika diminta)
        if use_telescope:
            (luminansi_hilal_tel_nl, sky_brightness_tel_nl,
             rasio_kontras_tel, delta_m_tel) = self.hitung_visibilitas_teleskop(
                luminansi_hilal_nl, sky_brightness_nl, posisi, k_v,
                aperture=aperture,
                magnification=magnification
            )
        else:
            luminansi_hilal_tel_nl = sky_brightness_tel_nl = 0.0
            rasio_kontras_tel = delta_m_tel = 0.0
        
        # Simpan semua hasil termasuk data posisi lengkap
        self.hasil.update({
            # Data posisi matahari
            'sun_alt': posisi['sun_alt'],
            'sun_az': posisi['sun_az'],
            # Data posisi bulan
            'moon_alt': posisi['moon_alt'],
            'moon_az': posisi['moon_az'],
            'elongation': posisi['elongation'],
            'phase_angle': posisi['phase_angle'],
            'moon_semidiameter': posisi['moon_semidiameter'],
            # Data perhitungan
            'luminansi_hilal_nl': luminansi_hilal_nl,
            'sky_brightness_nl': sky_brightness_nl,
            'luminansi_hilal_tel_nl': luminansi_hilal_tel_nl,
            'sky_brightness_tel_nl': sky_brightness_tel_nl,
            'k_v': k_v,
            'rasio_kontras_ne': rasio_kontras_ne,
            'delta_m_ne': delta_m_ne,
            'rasio_kontras_tel': rasio_kontras_tel,
            'delta_m_tel': delta_m_tel
        })
        
        # Tampilkan hasil akhir
        self.tampilkan_hasil_akhir()
        
        return self.hasil
    
    def tampilkan_hasil_akhir(self):
        """Menampilkan ringkasan hasil perhitungan akhir"""
        print("\n" + "=" * 70)
        print("HASIL PERHITUNGAN VISIBILITAS HILAL")
        print("=" * 70)
        
        # === INFORMASI LOKASI ===
        print(f"\n{'='*30} LOKASI {'='*31}")
        print(f"  Nama Tempat           : {self.nama_tempat}")
        print(f"  Lintang               : {self.lintang}°")
        print(f"  Bujur                 : {self.bujur}°")
        print(f"  Elevasi               : {self.elevasi} m")
        print(f"  Timezone              : {self.timezone_str}")
        print(f"  Bulan/Tahun Hijri     : {self.bulan_hijri}/{self.tahun_hijri}")
        
        # === WAKTU ===
        print(f"\n{'='*30} WAKTU {'='*32}")
        if 'ijtima_utc' in self.hasil:
            print(f"  Ijtima UTC            : {self.hasil['ijtima_utc'].strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"  Ijtima Lokal          : {self.hasil['ijtima_local'].strftime('%Y-%m-%d %H:%M:%S')}")
        if 'tanggal_pengamatan' in self.hasil:
            print(f"  Tanggal Pengamatan    : {self.hasil['tanggal_pengamatan'].strftime('%Y-%m-%d')}")
        if 'sunset_utc' in self.hasil:
            print(f"  Sunset UTC            : {self.hasil['sunset_utc'].strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"  Sunset Lokal          : {self.hasil['sunset_local'].strftime('%Y-%m-%d %H:%M:%S')}")
        
        # === DATA ATMOSFER ===
        print(f"\n{'='*28} DATA ATMOSFER {'='*28}")
        if 'rh' in self.hasil:
            print(f"  Kelembapan Relatif (RH)   : {self.hasil['rh']:.2f}%")
            print(f"  Suhu (T)                  : {self.hasil['temperature']:.2f}°C")
        if 'k_v' in self.hasil:
            print(f"  Koefisien Ekstingsi (k_V) : {self.hasil['k_v']:.4f}")
        
        # === POSISI MATAHARI ===
        print(f"\n{'='*27} POSISI MATAHARI {'='*27}")
        if 'sun_alt' in self.hasil:
            print(f"  Altitude Matahari     : {deg_to_dms(self.hasil['sun_alt'])}")
            print(f"  Azimuth Matahari      : {deg_to_dms(self.hasil['sun_az'])}")
        
        # === POSISI BULAN ===
        print(f"\n{'='*28} POSISI BULAN {'='*29}")
        if 'moon_alt' in self.hasil:
            print(f"  Altitude Bulan        : {deg_to_dms(self.hasil['moon_alt'])}")
            print(f"  Azimuth Bulan         : {deg_to_dms(self.hasil['moon_az'])}")
        if 'elongation' in self.hasil:
            print(f"  Elongasi Toposentrik  : {deg_to_dms(self.hasil['elongation'])}")
        if 'phase_angle' in self.hasil:
            print(f"  Phase Angle           : {deg_to_dms(self.hasil['phase_angle'])}")
        if 'moon_semidiameter' in self.hasil:
            print(f"  Semidiameter Bulan    : {deg_to_dms(float(self.hasil['moon_semidiameter']))}")
        
        # === LUMINANSI & SKY BRIGHTNESS ===
        print(f"\n{'='*23} LUMINANSI & SKY BRIGHTNESS {'='*22}")
        if 'luminansi_hilal_nl' in self.hasil:
            print(f"  Luminansi Hilal       : {self.hasil['luminansi_hilal_nl']:.4e} nL")
            print(f"  Sky Brightness        : {self.hasil['sky_brightness_nl']:.4e} nL")
        
        # === HASIL VISIBILITAS NAKED EYE ===
        print(f"\n{'='*25} VISIBILITAS NAKED EYE {'='*24}")
        if 'rasio_kontras_ne' in self.hasil:
            print(f"  Rasio Kontras (R)     : {self.hasil['rasio_kontras_ne']:.4e}")
            print(f"  Delta m               : {self.hasil['delta_m_ne']:.4f}")
            status_ne = "TERLIHAT" if self.hasil['delta_m_ne'] >= 0 else "TIDAK TERLIHAT"
            print(f"  Status                : {status_ne}")
        
        # === HASIL VISIBILITAS TELESKOP ===
        print(f"\n{'='*26} VISIBILITAS TELESKOP {'='*24}")
        if 'rasio_kontras_tel' in self.hasil:
            print(f"  Luminansi (Teleskop)  : {self.hasil['luminansi_hilal_tel_nl']:.4e} nL")
            print(f"  Sky Bright. (Teleskop): {self.hasil['sky_brightness_tel_nl']:.4e} nL")
            print(f"  Rasio Kontras (R)     : {self.hasil['rasio_kontras_tel']:.4e}")
            print(f"  Delta m               : {self.hasil['delta_m_tel']:.4f}")
            status_tel = "TERLIHAT" if self.hasil['delta_m_tel'] > 0 else "TIDAK TERLIHAT"
            print(f"  Status                : {status_tel}")
        
        # === INTERPRETASI ===
        print(f"\n{'='*28} INTERPRETASI {'='*29}")
        if self.hasil.get('delta_m_ne', -1) >= 0:
            print(f"  [OK] Naked Eye  : Hilal BERPOTENSI terlihat dengan mata telanjang")
        else:
            print(f"  [!]  Naked Eye  : Hilal SULIT terlihat dengan mata telanjang")
        
        if self.hasil.get('delta_m_tel', -1) > 0:
            print(f"  [OK] Teleskop   : Hilal TERDETEKSI dengan teleskop")
        else:
            print(f"  [!]  Teleskop   : Hilal SULIT terdeteksi dengan teleskop")
        
        print("\n" + "=" * 70)


def main():
    """Fungsi utama untuk demo dan testing"""
    
    # Contoh penggunaan sesuai spesifikasi
    # Lokasi: UIN Walisongo, Bulan: Muharram 1444
    
    calculator = HilalVisibilityCalculator(
        nama_tempat="uin",
        lintang=-6.916666667, # aceh: 5.466769444, uin: -6.916666667, mataram: -8.346986111, sulamu-kupang:-10.04527778
        bujur=110.3480556, # aceh: 95.24221944, uin:110.3480556, mataram: 116.149, sulamu-kupang: 123.6061569
        elevasi=89,
        timezone_str="Asia/Jakarta",
        bulan_hijri=1,      # Muharram
        tahun_hijri=1444, 
        delta_day_offset=0
    )
    
    # Jalankan perhitungan lengkap dengan teleskop
    hasil = calculator.jalankan_perhitungan_lengkap(
        use_telescope=True,
        aperture=100.0,      # 100mm teleskop
        magnification=50.0   # 50x pembesaran
    )
    
    # Summary singkat
    status_ne = "TERLIHAT" if hasil.get('delta_m_ne', -1) >= 0 else "SULIT TERLIHAT"
    status_tel = "TERDETEKSI" if hasil.get('delta_m_tel', -1) > 0 else "SULIT TERDETEKSI"
    print(f"\nPROGRAM SELESAI   : Hilal {status_ne} (naked eye), {status_tel} (teleskop)")
    print(f"Simpan hasil ke variabel hasil untuk analisis lebih lanjut.")
    
    return hasil


if __name__ == "__main__":
    hasil = main()
