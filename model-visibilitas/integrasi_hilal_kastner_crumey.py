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

1. IMPLEMENTASI MODEL CRUMEY 2014 UNTUK TELESKOP
   -----------------------------------------------
   Peningkatan: Mengganti model Schaefer dengan model Crumey untuk visibilitas teleskop:

   Model Baru (Crumey 2014):
   - Berdasarkan threshold contrast untuk extended source
   - Lebih realistis: TIDAK ada "boost" artefaktual M^2
   - Output: Visibility margin (log10(C/C_thr)), >0 = detectable
   - Fungsi: visibility_margin() dari crumey_telescope_correction.py

   Keuntungan:
   - Menghilangkan artefak "boost M^2" yang menyebabkan delta_m tidak wajar
   - Model berdasarkan human vision threshold (lebih ilmiah)
   - Single source of truth untuk luas sabit (crescent_area_arcmin2)

   Perubahan:
   - Import: from crumey_telescope_correction import visibility_margin, TelescopeSpec
   - Fungsi hitung_visibilitas_teleskop() sekarang menggunakan model Crumey
   - Output delta_m_tel sekarang adalah margin_tel (Crumey visibility margin)

   Catatan:
   - Naked eye tetap menggunakan model Schaefer (delta m)
   - Untuk extended source, luminansi teleskop SAMA dengan naked eye (tidak ada boost)
   - Margin > 0 berarti contrast objek melebihi threshold deteksi human eye

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
# from limit_teleskop import extended_surface_correction_factor  # REPLACED with Crumey
from crumey_telescope_correction import visibility_margin, TelescopeSpec, crescent_area_arcmin2
from meteo_power_api import ObservingLocation, get_rh_t_at_time, PowerAPIError

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
        print("Menghitung waktu ijtima/konjungsi...")
        ijtima_utc = newmoon_hijri_month_utc(self.tahun_hijri, self.bulan_hijri)
        ijtima_local = convert_utc_to_localtime(self.timezone_str, utc_datetime=ijtima_utc)
        
        self.hasil['ijtima_utc'] = ijtima_utc
        self.hasil['ijtima_local'] = ijtima_local
        
        print(f"  Ijtima UTC  : {ijtima_utc.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Ijtima Lokal: {ijtima_local.strftime('%Y-%m-%d %H:%M:%S')}")
        
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
        print("\nMenentukan tanggal pengamatan...")

        # Konversi ijtima ke waktu lokal
        ijtima_local = convert_utc_to_localtime(self.timezone_str, utc_datetime=ijtima_utc)

        # Aturan sederhana berdasarkan jam ijtima lokal
        if ijtima_local.hour < 12:
            # Ijtima terjadi sebelum jam 12:00 (tengah malam sampai sebelum siang)
            # Pengamatan dilakukan pada hari yang sama (tanggal ijtima lokal)
            tanggal_pengamatan = ijtima_local.date()
            print(f"  Ijtima lokal sebelum jam 12:00")
            print(f"  Ijtima Lokal: {ijtima_local.strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"  Pengamatan  : {tanggal_pengamatan.strftime('%Y-%m-%d')} (hari ijtima lokal, sore hari)")
        else:
            # Ijtima terjadi setelah jam 12:00 (siang sampai tengah malam)
            # Hitung sunset pada hari ijtima lokal (gunakan sunset geometris untuk perbandingan)
            _, sunset_local_ijtima_geom = sunrise_sunset_geometric_local(
                self.location,
                self.timezone_str,
                year=ijtima_utc.year,
                month=ijtima_utc.month,
                day=ijtima_utc.day
            )

            # Bandingkan ijtima lokal dengan sunset lokal
            if ijtima_local < sunset_local_ijtima_geom:
                # Ijtima sebelum sunset sore, amati hari berikutnya
                tanggal_pengamatan = ijtima_local.date() + timedelta(days=1)
                print(f"  Ijtima lokal setelah jam 12:00 dan sebelum sunset sore")
                print(f"  Ijtima Lokal : {ijtima_local.strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"  Sunset Lokal : {sunset_local_ijtima_geom.strftime('%Y-%m-%d %H:%M:%S')} (geometris)")
                print(f"  Pengamatan   : {tanggal_pengamatan.strftime('%Y-%m-%d')} (hari berikutnya)")
            else:
                # Ijtima setelah sunset sore, amati hari yang sama
                tanggal_pengamatan = ijtima_local.date()
                print(f"  Ijtima lokal setelah jam 12:00 dan setelah sunset sore")
                print(f"  Ijtima Lokal : {ijtima_local.strftime('%Y-%m-%d %H:%M:%S')}")
                print(f"  Sunset Lokal : {sunset_local_ijtima_geom.strftime('%Y-%m-%d %H:%M:%S')} (geometris)")
                print(f"  Pengamatan   : {tanggal_pengamatan.strftime('%Y-%m-%d')} (hari yang sama)")

        # Terapkan delta day offset
        tanggal_pengamatan += timedelta(days=self.delta_day_offset)

        # LANGKAH 1: Hitung sunset GEOMETRIS (tanpa koreksi refraksi)
        print(f"\n  Menghitung sunset geometris (tanpa koreksi refraksi)...")
        _, sunset_geom_local = sunrise_sunset_geometric_local(
            self.location,
            self.timezone_str,
            year=tanggal_pengamatan.year,
            month=tanggal_pengamatan.month,
            day=tanggal_pengamatan.day
        )
        sunset_geom_utc = convert_localtime_to_utc(self.timezone_str, local_datetime=sunset_geom_local)
        print(f"  Sunset Geometris UTC  : {sunset_geom_utc.strftime('%H:%M:%S')}")
        print(f"  Sunset Geometris Lokal: {sunset_geom_local.strftime('%H:%M:%S')}")

        # LANGKAH 2: Ambil data atmosfer (RH dan T) pada waktu sunset geometris
        print(f"\n  Mengambil data atmosfer dari NASA POWER API...")
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
            print(f"  Kelembapan Relatif (RH): {rh:.2f}%")
            print(f"  Suhu (T): {temperature:.2f} derajat Celsius")
        except (PowerAPIError, Exception) as e:
            print(f"  Warning: Gagal mengambil data dari NASA POWER API ({e})")
            print(f"  Menggunakan nilai default: RH=80%, T=25°C")
            rh = 80.0
            temperature = 25.0

        # LANGKAH 3: Hitung sunset APPARENT dengan koreksi refraksi menggunakan T dari API
        print(f"\n  Menghitung sunset apparent dengan koreksi refraksi atmosfer...")
        print(f"  Menggunakan suhu: {temperature:.2f}°C dari API")
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
        print(f"  Sunset Apparent UTC  : {sunset_utc.strftime('%H:%M:%S')}")
        print(f"  Sunset Apparent Lokal: {sunset_local.strftime('%H:%M:%S')}")

        # Simpan hasil
        self.hasil['sunset_utc'] = sunset_utc
        self.hasil['sunset_local'] = sunset_local
        self.hasil['sunset_geom_utc'] = sunset_geom_utc
        self.hasil['sunset_geom_local'] = sunset_geom_local
        self.hasil['tanggal_pengamatan'] = datetime.combine(tanggal_pengamatan, datetime.min.time())
        self.hasil['rh'] = rh
        self.hasil['temperature'] = temperature

        print(f"\n  Tanggal Pengamatan: {tanggal_pengamatan.strftime('%Y-%m-%d')}")
        print(f"  Sunset Final (dengan refraksi): {sunset_local.strftime('%H:%M:%S')}")

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
        print("\nMengambil data atmosfer dari NASA POWER API...")
        
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
            # Jika naive, anggap sudah dalam UTC dan beri tzinfo
            sunset_utc = sunset_utc.replace(tzinfo=timezone.utc)
        else:
            # Jika bertimezone tapi bukan UTC, konversi ke UTC
            sunset_utc = sunset_utc.astimezone(timezone.utc)

        try:
            rh, temperature = get_rh_t_at_time(loc, sunset_utc)
            print(f"  Kelembapan Relatif (RH): {rh:.2f}%")
            print(f"  Suhu (T): {temperature:.2f} derajat Celsius")
            
            self.hasil['rh'] = rh
            self.hasil['temperature'] = temperature
            
            return rh, temperature
        except PowerAPIError as e:
            print(f"  Error: Gagal mengambil data dari NASA POWER API")
            print(f"  Menggunakan nilai default: RH=80%, T=25°C")
            print(f"  Detail error: {e}")
            
            # Nilai default jika API gagal
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
        print("\nMenghitung posisi matahari dan bulan saat sunset...")
        
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
        
        print(f"  Sun Altitude : {sun_alt} derajat")
        print(f"  Sun Azimuth  : {sun_az} derajat")
        print(f"  Moon Altitude: {moon_alt} derajat")
        print(f"  Moon Azimuth : {moon_az} derajat")
        print(f"  Elongasi     : {elongation} derajat")
        print(f"  Phase Angle  : {phase_angle} derajat")
        print(f"  Moon Semidiam: {float(moon_semidiameter)} derajat")
        
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
        print("\nMenghitung sky brightness (model Schaefer)...")
        
        # Hitung selisih azimuth sun-moon
        azisun = abs(posisi['sun_az'] - posisi['moon_az'])
        
        # Gunakan nilai sun_alt apa adanya untuk keakuratan maksimal
        # Tidak ada pembulatan atau batasan nilai
        sun_alt_safe = posisi['sun_alt']
        
        try:
            # Panggil fungsi visual_limit
            result = visual_limit(
                month=self.hasil['tanggal_pengamatan'].month,
                year=self.hasil['tanggal_pengamatan'].year,
                phase_angle=self.PHASE_MOON,  # Statis: 0 (new moon)
                altmoon=self.ALTMOON_STATIS,  # Statis: -90 derajat
                azimoon=self.AZIMOON_STATIS,  # Statis: 180 derajat
                altsun=sun_alt_safe,          # Otomatis (dengan batas aman)
                azisun=azisun,                 # Otomatis: selisih azimuth sun-moon
                humidity=rh,
                temperature=temperature,
                latitude=self.lintang,
                altitude=self.elevasi,
                snellen_ratio=self.SNELLEN_RATIO,  # Statis: 1
                altstar=max(posisi['moon_alt'], 0.0),  # Otomatis: altitude hilal (min 0)
                goal_magnitude=99
            )
            
            sky_brightness_nl = float(result.get("sky_brightness", 0.0))
            
            # Jika sky brightness negatif atau sangat kecil, gunakan nilai minimum
            if sky_brightness_nl <= 0:
                print(f"  Warning: Sky brightness negatif ({sky_brightness_nl:.6e})")
                print(f"  Menggunakan nilai minimum: 1.0 nL")
                sky_brightness_nl = 1.0
            
        except (ValueError, ZeroDivisionError) as e:
            print(f"  Error dalam perhitungan sky brightness: {e}")
            print(f"  Menggunakan nilai default: 1000.0 nL")
            sky_brightness_nl = 1000.0
        
        # Ambil koefisien ekstingsi k_V dari array K (index 2 untuk V-band)
        K_list = result.get("K", []) if 'result' in locals() else []
        if len(K_list) > 2:
            k_v = float(K_list[2])
        else:
            k_v = 0.3  # Default value
        
        print(f"  Sky Brightness (nL): {sky_brightness_nl}")
        print(f"  Koefisien Ekstingsi (k_V): {k_v}")
        
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
        print("\nMenghitung luminansi hilal (model Kastner)...")
        
        # Zenith distance = 90 - altitude
        zenith_distance = 90.0 - posisi['moon_alt']
        
        # Hitung luminansi hilal
        luminansi_hilal_nl = hitung_luminansi_kastner(
            alpha=posisi['phase_angle'],
            elongation=posisi['elongation'],
            r=float(posisi['moon_semidiameter']),
            z=zenith_distance,
            k=k_v
        )
        
        print(f"  Luminansi Hilal (nL): {luminansi_hilal_nl}")
        
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
        print("\nMenghitung visibilitas naked eye...")
        
        if sky_brightness_nl <= 0.0:
            raise ValueError("Sky brightness harus positif")
        
        # Rasio kontras
        rasio_kontras = luminansi_hilal_nl / sky_brightness_nl
        
        # Delta m
        if rasio_kontras > 0.0:
            delta_m = 2.5 * math.log10(rasio_kontras)
        else:
            delta_m = float('-inf')
        
        print(f"  Rasio Kontras (R): {rasio_kontras}")
        print(f"  Delta m: {delta_m}")
        
        return rasio_kontras, delta_m
    
    def hitung_visibilitas_teleskop(self,
                                     luminansi_hilal_nl: float,
                                     sky_brightness_nl: float,
                                     posisi: Dict[str, float],
                                     aperture: float = 100.0,
                                     magnification: float = 50.0,
                                     transmission: float = 0.95,
                                     n_surfaces: int = 6,
                                     central_obstruction: float = 0.0,
                                     observer_age: float = 22.0,
                                     field_factor: float = 2.4) -> Tuple[float, float, float, float]:
        """
        Menghitung visibilitas untuk pengamatan dengan teleskop menggunakan model Crumey (2014).

        Model Crumey berdasarkan threshold contrast untuk extended source, lebih realistis
        daripada pendekatan Schaefer yang memberikan "boost" artefaktual M^2.

        Parameters:
        -----------
        luminansi_hilal_nl : float
            Luminansi hilal dalam nL
        sky_brightness_nl : float
            Sky brightness dalam nL
        posisi : Dict[str, float]
            Dictionary posisi matahari dan bulan (untuk elongation dan moon_semidiameter)
        aperture : float
            Diameter aperture teleskop (mm)
        magnification : float
            Pembesaran teleskop
        transmission : float
            Transmisi per permukaan optik
        n_surfaces : int
            Jumlah permukaan optik
        central_obstruction : float
            Diameter obstruksi pusat (mm)
        observer_age : float
            Usia pengamat (tahun)
        field_factor : float
            Faktor konservatif untuk kondisi lapangan (default: 2.4)

        Returns:
        --------
        luminansi_hilal_tel_nl : float
            Luminansi hilal teleskop dalam nL (sama dengan input, no boost)
        sky_brightness_tel_nl : float
            Sky brightness teleskop dalam nL (sama dengan input untuk extended source)
        rasio_kontras_tel : float
            Rasio kontras teleskop (sama dengan naked eye untuk extended source)
        margin_tel : float
            Margin visibilitas teleskop (log10(C/C_thr)), >0 = detectable
        """
        print("\nMenghitung visibilitas teleskop (Model Crumey 2014)...")

        # Hitung transmisi per permukaan
        if n_surfaces > 0 and transmission > 0:
            t_per_surface = transmission ** (1.0 / n_surfaces)
        else:
            t_per_surface = transmission

        # Hitung luas sabit bulan dalam arcmin²
        elongation = posisi['elongation']
        r_deg = float(posisi['moon_semidiameter'])
        A_arcmin2 = crescent_area_arcmin2(elongation, r_deg)

        print(f"  Elongasi: {elongation:.4f} derajat")
        print(f"  Semidiameter: {r_deg:.6f} derajat")
        print(f"  Luas Sabit (A): {A_arcmin2:.4f} arcmin²")

        # Buat TelescopeSpec untuk Crumey
        spec = TelescopeSpec(
            aperture_mm=aperture,
            magnification=magnification,
            obstruction_mm=central_obstruction,
            transmission_per_surface=t_per_surface,
            n_surfaces=n_surfaces,
            observer_age=observer_age,
            field_factor=field_factor,
            phi=1.0
        )

        # Hitung visibility margin dengan model Crumey
        out = visibility_margin(luminansi_hilal_nl, sky_brightness_nl, A_arcmin2,
                               spec=spec, use_abs_contrast=True)

        # Untuk extended source, contrast Weber (L-B)/B TIDAK berubah dengan magnifikasi
        # Jadi luminansi_hilal_tel dan sky_brightness_tel sama dengan input
        luminansi_hilal_tel_nl = luminansi_hilal_nl
        sky_brightness_tel_nl = sky_brightness_nl

        # Rasio kontras sama (untuk extended source)
        if sky_brightness_tel_nl > 0.0:
            rasio_kontras_tel = luminansi_hilal_tel_nl / sky_brightness_tel_nl
        else:
            rasio_kontras_tel = float('nan')

        # Margin visibilitas teleskop dari Crumey
        margin_tel = out['margin_tel']
        C_thr_tel = out['C_thr_tel']

        print(f"  Luminansi Hilal (nL): {luminansi_hilal_tel_nl:.2f}")
        print(f"  Sky Brightness (nL): {sky_brightness_tel_nl:.2f}")
        print(f"  Rasio Kontras (R): {rasio_kontras_tel:.2f}")
        print(f"  Threshold Contrast (C_thr): {C_thr_tel:.6f}")
        print(f"  Margin Teleskop: {margin_tel:.4f}")
        print(f"  Status: {'DETECTABLE' if margin_tel > 0 else 'NOT DETECTABLE'}")

        # Simpan data diagnostic dari Crumey
        self.hasil['crumey_diagnostic'] = {
            'A_arcmin2': A_arcmin2,
            'Ba_cd_m2': out.get('Ba_cd_m2', 0),
            'C_thr_tel': C_thr_tel,
            'exit_pupil_mm': out.get('exit_pupil_mm', 0),
            'throughput_T': out.get('throughput_T', 0),
        }

        # Return dengan margin_tel sebagai pengganti delta_m_tel
        return luminansi_hilal_tel_nl, sky_brightness_tel_nl, rasio_kontras_tel, margin_tel
    
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
        print("=" * 70)
        print("PROGRAM INTEGRASI MODEL VISIBILITAS HILAL KASTNER MODIFIKASI")
        print("=" * 70)
        print(f"\nLokasi: {self.nama_tempat}")
        print(f"Koordinat: Lintang={self.lintang} derajat, Bujur={self.bujur} derajat, Elevasi={self.elevasi} m")
        print(f"Timezone: {self.timezone_str}")
        print(f"Bulan Hijri: {self.bulan_hijri} Tahun {self.tahun_hijri}")
        print(f"Delta Day Offset: {self.delta_day_offset}")

        # Langkah 1: Hitung ijtima
        ijtima_utc, ijtima_local = self.hitung_ijtima()

        # Langkah 2: Tentukan tanggal pengamatan dan hitung sunset dengan koreksi refraksi
        # Fungsi ini sudah mengambil data atmosfer dan menghitung sunset dengan alur yang benar
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
             rasio_kontras_tel, margin_tel) = self.hitung_visibilitas_teleskop(
                luminansi_hilal_nl, sky_brightness_nl, posisi,
                aperture=aperture,
                magnification=magnification
            )
            # Untuk backward compatibility, gunakan margin_tel sebagai delta_m_tel
            delta_m_tel = margin_tel
        else:
            luminansi_hilal_tel_nl = sky_brightness_tel_nl = 0.0
            rasio_kontras_tel = delta_m_tel = 0.0
        
        # Simpan semua hasil
        self.hasil.update({
            'moon_alt': posisi['moon_alt'],
            'elongation': posisi['elongation'],
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
        
        print(f"\nInformasi Dasar:")
        print(f"  Nama Tempat          : {self.nama_tempat}")
        print(f"  Koordinat            : Lintang={self.lintang} derajat, Bujur={self.bujur} derajat, Elevasi={self.elevasi} m")
        print(f"  Time Zone            : {self.timezone_str}")
        
        if 'ijtima_utc' in self.hasil:
            print(f"\nWaktu Ijtima:")
            print(f"  Ijtima UTC           : {self.hasil['ijtima_utc'].strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"  Ijtima Local         : {self.hasil['ijtima_local'].strftime('%Y-%m-%d %H:%M:%S')}")
        
        if 'sunset_utc' in self.hasil:
            print(f"\nWaktu Sunset:")
            print(f"  Sunset UTC           : {self.hasil['sunset_utc'].strftime('%Y-%m-%d %H:%M:%S')}")
            print(f"  Sunset Local         : {self.hasil['sunset_local'].strftime('%Y-%m-%d %H:%M:%S')}")
        
        print(f"\nDelta Day Offset       : {self.delta_day_offset}")
        
        if 'moon_alt' in self.hasil:
            print(f"\nData Posisi Bulan:")
            print(f"  Moon Altitude        : {self.hasil['moon_alt']} derajat")
            print(f"  Elongasi Toposentrik : {self.hasil['elongation']} derajat")
        
        if 'rh' in self.hasil:
            print(f"\nData Atmosfer:")
            print(f"  Relative Humidity (RH): {self.hasil['rh']}%")
            print(f"  Temperature (T)      : {self.hasil['temperature']} derajat Celsius")
        
        if 'luminansi_hilal_nl' in self.hasil:
            print(f"\nHasil Perhitungan:")
            print(f"  Luminansi Hilal (Naked Eye)     : {self.hasil['luminansi_hilal_nl']} nL")
            print(f"  Sky Brightness (Naked Eye)      : {self.hasil['sky_brightness_nl']} nL")
            print(f"  Luminansi Hilal (Teleskop)      : {self.hasil['luminansi_hilal_tel_nl']} nL")
            print(f"  Sky Brightness (Teleskop)       : {self.hasil['sky_brightness_tel_nl']} nL")
            print(f"  Koefisien Ekstingsi (k_V)       : {self.hasil['k_v']}")
            print(f"  Rasio Kontras R (Naked Eye)     : {self.hasil['rasio_kontras_ne']}")
            print(f"  Rasio Kontras R (Teleskop)      : {self.hasil['rasio_kontras_tel']}")
            print(f"  Delta m (Naked Eye)             : {self.hasil['delta_m_ne']}")
            print(f"  Margin Teleskop (Crumey)        : {self.hasil['delta_m_tel']}")

        # Tampilkan diagnostic Crumey jika tersedia
        if 'crumey_diagnostic' in self.hasil:
            diag = self.hasil['crumey_diagnostic']
            print(f"\nDiagnostic Crumey:")
            print(f"  Luas Sabit (A_arcmin2)          : {diag.get('A_arcmin2', 0):.4f} arcmin²")
            print(f"  Exit Pupil                      : {diag.get('exit_pupil_mm', 0):.2f} mm")
            print(f"  Throughput                      : {diag.get('throughput_T', 0):.4f}")
            print(f"  Threshold Contrast (C_thr)      : {diag.get('C_thr_tel', 0):.6f}")

        # Interpretasi hasil
        print(f"\nInterpretasi:")
        if self.hasil['delta_m_ne'] < 0:
            print(f"  [!] Delta m (Naked Eye) < 0  -> Hilal SULIT untuk dilihat dengan mata telanjang")
        else:
            print(f"  [OK] Delta m (Naked Eye) >= 0  -> Hilal mungkin terlihat dengan mata telanjang")

        # Untuk teleskop, gunakan margin dari Crumey (>0 = detectable)
        if self.hasil['delta_m_tel'] < 0:
            print(f"  [!] Margin Teleskop < 0  -> Hilal SULIT untuk dilihat dengan teleskop")
        else:
            print(f"  [OK] Margin Teleskop > 0  -> Hilal TERDETEKSI dengan teleskop (Model Crumey)")

        # Catatan tentang model
        print(f"\nCatatan:")
        print(f"  - Naked eye menggunakan model Schaefer (delta m)")
        print(f"  - Teleskop menggunakan model Crumey 2014 (visibility margin)")
        print(f"  - Margin > 0 berarti contrast objek melebihi threshold deteksi")
        
        print("\n" + "=" * 70)


def main():
    """Fungsi utama untuk demo dan testing"""
    
    # Contoh penggunaan sesuai spesifikasi
    # Lokasi: UIN Walisongo
    # Bulan: Muharram 1444
    
    calculator = HilalVisibilityCalculator(
        nama_tempat="UIN Walisongo Semarang",
        lintang=-6.9666,
        bujur=110.45,
        elevasi=89.0,
        timezone_str="Asia/Jakarta",
        bulan_hijri=7,      # Muharram
        tahun_hijri=1447,
        delta_day_offset=0
    )
    
    # Jalankan perhitungan lengkap dengan teleskop
    hasil = calculator.jalankan_perhitungan_lengkap(
        use_telescope=True,
        aperture=100.0,      # 100mm teleskop
        magnification=50.0   # 50x pembesaran
    )
    
    print("\n" + "="*70)
    print("PROGRAM SELESAI")
    print("="*70)
    
    return hasil


if __name__ == "__main__":
    hasil = main()
    print(f"\nSimpan hasil ke variabel hasil untuk analisis lebih lanjut.")
