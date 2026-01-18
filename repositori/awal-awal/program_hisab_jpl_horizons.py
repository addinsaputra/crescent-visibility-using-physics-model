"""
Program Hisab Hilal menggunakan data dari JPL Horizons
======================================================

File ini menyediakan serangkaian fungsi untuk mengambil dan
mengolah data ephemeris dari layanan **JPL Horizons** sebagai
dokumen Python. Program ini dirancang untuk membantu proses
``hisab`` awal bulan Hijriah dengan model Schaefer serta
modifikasi model Kastner. Modul ini tidak melakukan perhitungan
penetapan hilal secara langsung, melainkan menyediakan utilitas
untuk:

* Menemukan waktu konjungsi Matahari–Bulan (ijtima') secara
  geosentris dengan mencari minimum sudut elongasi Bulan dari
  Matahari. Sudut elongasi didefinisikan oleh Horizons sebagai
  **S‑O‑T** (Sun–Observer–Target)【482484488989631†L733-L736】.
* Mengambil waktu matahari terbenam pada lokasi tertentu untuk
  tanggal tertentu melalui opsi **R_T_S_ONLY=YES** di API Horizons.
* Mengambil azimut dan ketinggian (elevasi) Matahari dan Bulan
  dari lokasi topocentrik. Parameter *quantities* bernilai **4**
  mengembalikan kolom `AZ` (Azimuth) dan `EL` (Elevation)
  【482484488989631†L580-L589】.
* Mengambil sudut elongasi (S‑O‑T), sudut fase (phi), diameter
  sudut Bulan, fraksi iluminasi, dan jarak Bulan dari pengamat
  (delta)【482484488989631†L733-L744】. Kolom‐kolom ini diperlukan
  sebagai input model Kastner.

Untuk menjalankan skrip ini diperlukan paket `requests` (sudah
terinstal pada lingkungan ini) serta modul `zoneinfo` dari
Python 3.9+ untuk manajemen zona waktu. Jika Anda ingin
menggunakan `astroquery` sebagai pembungkus API Horizons, modul
tersebut harus diinstal secara manual (misalnya `pip install
astroquery`). Skrip ini menggunakan pendekatan HTTP murni agar
tidak bergantung pada modul eksternal.

Catatan: Contoh kode ini menulis fungsi yang memanggil API
Horizons. Agar berjalan di lingkungan Anda, koneksi internet
dibutuhkan. Jika server `ssd.jpl.nasa.gov` sedang mengalami
perubahan sertifikat, pertimbangkan untuk menambahkan parameter
`verify=False` pada pemanggilan `requests.get` (walaupun hal ini
tidak direkomendasikan di luar debugging karena menonaktifkan
verifikasi TLS).

Referensi API: Parameter `quantities` pada `ephemerides` memberi
kesempatan untuk memilih subset kolom sesuai definisi “Observer
Table Quantities” di manual Horizons【997971877924426†L294-L303】.

"""

from __future__ import annotations

import json
import math
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

import requests
from zoneinfo import ZoneInfo

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class HorizonsClient:
    """
    Klien sederhana untuk berinteraksi dengan API JPL Horizons
    menggunakan permintaan HTTP. Objek ini menampung informasi
    lokasi pengamat (topocentrik) dan zona waktu lokal, serta
    menyediakan fungsi untuk mengambil ephemeris Matahari, Bulan,
    atau objek lain.

    Attributes
    ----------
    location : dict
        Kamus berisi `lon` (longitude derajat), `lat` (latitude derajat),
        dan `elevation` (ketinggian km). Nilai longitude harus
        mengikuti konvensi Horizons: derajat timur positif untuk
        Bumi, sehingga lokasi 110° BT direpresentasikan sebagai
        `110.0`.
    timezone : ZoneInfo
        Zona waktu lokal pengamat. Digunakan untuk konversi waktu.
    base_url : str
        Alamat dasar API Horizons. Secara default
        ``https://ssd.jpl.nasa.gov/api/horizons.api``.
    """

    def __init__(self, location: Dict[str, float], timezone: str = "UTC") -> None:
        self.location = location
        self.timezone = ZoneInfo(timezone)
        self.base_url = "https://ssd.jpl.nasa.gov/api/horizons.api"

    def _make_request(self, params: Dict[str, str]) -> Dict:
        """
        Membuat permintaan GET ke API Horizons.

        Parameter `params` berisi semua parameter query sesuai
        dokumentasi Horizons. Fungsi mengembalikan dict hasil JSON
        apabila parameter `format=json` diberikan; jika tidak,
        mengembalikan dict dengan kunci ``text`` yang memuat string
        response mentah.
        """
        logger.debug("Memanggil API Horizons dengan parameter: %s", params)
        response = requests.get(self.base_url, params=params, timeout=60)
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            logger.error("HTTP Error: %s", e)
            logger.error("Response Text: %s", response.text)
            raise
        # Coba parse JSON; jika gagal, kembalikan sebagai teks.
        try:
            return response.json()
        except json.JSONDecodeError:
            return {"text": response.text}

    def _build_site_coord(self) -> str:
        """
        Membangun string `SITE_COORD` dari kamus lokasi. Horizons
        menerima format "lon,lat,elev" dengan satuan derajat dan km.
        """
        return f"{self.location['lon']},{self.location['lat']},{self.location['elevation']}"

    def get_sunrise_sunset(self, date: datetime) -> Tuple[datetime, datetime]:
        """
        Mengambil waktu terbit dan terbenam Matahari untuk tanggal
        tertentu pada lokasi topocentrik. Horizons menyediakan
        parameter ``R_T_S_ONLY=YES`` untuk mengembalikan daftar waktu
        terbit (rise), transit, dan terbenam (set).

        Parameters
        ----------
        date : datetime
            Tanggal (tanpa jam) dalam zona waktu lokal yang ingin
            ditentukan waktu terbit/terbenamnya. Waktu yang
            diberikan jamnya diabaikan; hanya tanggal digunakan.

        Returns
        -------
        tuple
            Pasangan ``(terbit_utc, terbenam_utc)`` dalam zona waktu
            UTC. Bila Horizons tidak mengembalikan hasil, keduanya
            bernilai ``None``.
        """
        # Format tanggal sebagai rentang 1 hari penuh dalam UTC.
        local_date = datetime(date.year, date.month, date.day, tzinfo=self.timezone)
        start_local = local_date
        stop_local = local_date + timedelta(days=1)
        params = {
            "format": "json",
            "COMMAND": "'10'",  # 10 = Sun
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "OBSERVER",
            "CENTER": "'coord@399'",  # topocentric Earth
            "COORD_TYPE": "'GEODETIC'",
            "SITE_COORD": f"'{self._build_site_coord()}'",
            "START_TIME": f"'{start_local.astimezone(ZoneInfo('UTC')).strftime('%Y-%m-%d %H:%M')}'",
            "STOP_TIME": f"'{stop_local.astimezone(ZoneInfo('UTC')).strftime('%Y-%m-%d %H:%M')}'",
            "STEP_SIZE": "'1 m'",
            "R_T_S_ONLY": "'YES'",
            "CSV_FORMAT": "'YES'",
        }
        result = self._make_request(params)
        rises, sets = None, None
        # API JSON menempatkan data ephemeris di bawah kunci 'result'
        if isinstance(result, dict) and 'result' in result:
            lines = result['result'].splitlines()
            in_data = False
            for line in lines:
                if '$$SOE' in line:
                    in_data = True
                    continue
                if '$$EOE' in line:
                    break
                if not in_data:
                    continue
                
                # Baris data format: YYYY-MM-DD HH:MM, [marker], Event, ...
                parts = line.strip().split(',')
                if len(parts) >= 3:
                    time_str = parts[0].strip()
                    # Event markers can be in col 1 or 2. e.g. "...,C,s,..."
                    p1 = parts[1].strip()
                    p2 = parts[2].strip() if len(parts) > 2 else ""
                    
                    event_code = None
                    # Prioritize finding 'R'/'S'/'r'/'s'
                    if p1.upper() in ['R', 'S']:
                        event_code = p1.upper()
                    elif p2.upper() in ['R', 'S']:
                        event_code = p2.upper()
                    
                    if not event_code:
                        continue

                    try:
                        # Format Horizons bisa "YYYY-Mon-DD HH:MM" -> need parsing flexibility
                        dt = datetime.strptime(time_str, '%Y-%b-%d %H:%M').replace(tzinfo=ZoneInfo('UTC'))
                        if event_code == 'R' and rises is None:
                            rises = dt
                        elif event_code == 'S' and sets is None:
                            sets = dt
                    except ValueError:
                        continue
        return rises, sets

    def query_ephemeris(self,
                        target_id: str,
                        times_utc: List[datetime],
                        quantities: str = '4,10,13,20,23,43') -> List[Dict[str, str]]:
        """
        Mengambil ephemeris untuk daftar waktu tertentu.

        Parameters
        ----------
        target_id : str
            Kode objek di Horizons. Misal `'301'` untuk Bulan, `'10'` untuk
            Matahari. Dapat juga diberikan sebagai nama ('Moon', 'Sun').
        times_utc : list of datetime
            Daftar waktu (UTC) yang akan dievaluasi.
        quantities : str, optional
            Daftar kode kuantitas Horizons yang ingin diambil. Default
            `'4,10,13,20,23,43'` mencakup azimut & elevasi, fraksi
            iluminasi, diameter sudut, jarak, elongasi, dan sudut fase
            true.

        Returns
        -------
        list of dict
            Daftar baris ephemeris. Masing‑masing baris adalah kamus
            dengan nama kolom sebagai kunci. Kolom `datetime_str` sudah
            dikonversi ke objek `datetime` dengan tzinfo UTC.
        """
        if not times_utc:
            return []
        tlist_str = ','.join([dt.strftime('%Y-%m-%d %H:%M') for dt in times_utc])
        params = {
            "format": "json",
            "COMMAND": f"'{target_id}'",
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "OBSERVER",
            "CENTER": "'coord@399'",
            "COORD_TYPE": "'GEODETIC'",
            "SITE_COORD": f"'{self._build_site_coord()}'",
            "TLIST": f"'{tlist_str}'",
            "TLIST_TYPE": "'CAL'",
            "TIME_TYPE": "'UT'",
            "CSV_FORMAT": "'YES'",
            "QUANTITIES": f"'{quantities}'",
            # Gunakan apparent refracted koordinat untuk visibilitas【482484488989631†L580-L589】
            "APPARENT": "'REFRACTED'",
        }
        result = self._make_request(params)
        rows: List[Dict[str, str]] = []
        if isinstance(result, dict) and 'result' in result:
            lines = result['result'].splitlines()
            header = None
            for line in lines:
                if not line.strip():
                    continue
                if header is None:
                    # Header baris biasanya mengandung nama kolom target like Azimuth, Elevation
                    # Horizons columns: Azi_(r-app), Elev_(r-app)
                    if ('Azi' in line or 'AZ' in line) and ('Elev' in line or 'EL' in line):
                        header = [col.strip() for col in line.split(',')]
                    continue
                parts = [p.strip() for p in line.split(',')]
                if header and len(parts) == len(header):
                    row = dict(zip(header, parts))
                    # Gunakan kolom pertama sebagai datetime
                    date_col = header[0]
                    if date_col in row:
                        dt_str = row[date_col].strip()
                        # Clean up " C" markers if present in date string? unique to some outputs?
                        # Usually simple split works.
                        # Handle potential seconds: 2026-May-11 10:30:00.000
                        try:
                            row['datetime_str'] = datetime.strptime(dt_str, '%Y-%b-%d %H:%M:%S.%f').replace(tzinfo=ZoneInfo('UTC'))
                        except ValueError:
                            try:
                                row['datetime_str'] = datetime.strptime(dt_str, '%Y-%b-%d %H:%M').replace(tzinfo=ZoneInfo('UTC'))
                            except ValueError:
                                pass
                    rows.append(row)
        return rows

    def find_geocentric_conjunction(self,
                                    start_utc: datetime,
                                    end_utc: datetime,
                                    step: timedelta = timedelta(minutes=10)) -> Optional[datetime]:
        """
        Mencari waktu ijtima' geosentris (minimum elongasi Bulan–Matahari) dalam
        interval yang diberikan. Fungsi ini mensample sudut elongasi
        (**S-O-T**) pada interval teratur lalu mencari nilai minimum.

        Parameters
        ----------
        start_utc, end_utc : datetime
            Batas awal dan akhir pencarian dalam UTC.
        step : timedelta
            Langkah sampling. Default setiap 10 menit. Anda dapat
            memperkecil langkah untuk presisi lebih tinggi.

        Returns
        -------
        datetime or None
            Waktu UTC saat sudut elongasi minimum ditemukan. Jika tidak
            ada data, mengembalikan ``None``.
        """
        # Gunakan pusat bumi sebagai pengamat: 500@399
        center = "'500@399'"
        # Helper untuk converting step to string (e.g. '10 m')
        step_minutes = int(step.total_seconds() // 60)
        
        params = {
            "format": "json",
            "COMMAND": "'301'",
            "MAKE_EPHEM": "YES",
            "EPHEM_TYPE": "OBSERVER",
            "CENTER": center,
            "START_TIME": f"'{start_utc.strftime('%Y-%m-%d %H:%M')}'",
            "STOP_TIME": f"'{end_utc.strftime('%Y-%m-%d %H:%M')}'",
            "STEP_SIZE": f"'{step_minutes} m'",
            "CSV_FORMAT": "'YES'",
            "QUANTITIES": "'23,31'",  # elongasi (23) dan bujur ekliptika (31)
            "APPARENT": "'AIRLESS'",
        }
        result = self._make_request(params)
        min_elon = 1e9
        min_time = None
        if isinstance(result, dict) and 'result' in result:
            lines = result['result'].splitlines()
            header = None
            for i, line in enumerate(lines):
                if i < 20: logger.info("RAW LINE %d: %s", i, line)
                if not line.strip():
                    continue
                # Skip preamble/metadata lines until we find the header
                if header is None:
                    # Header baris biasanya mengandung nama kolom target 'S-O-T'
                    if 'S-O-T' in line:
                         header = [col.strip() for col in line.split(',')]
                         # Hapus 'ERROR' jika ada di header karena kadang Horizons menyisipkannya
                         logger.debug("Found header: %s", header)
                    continue
                
                # Parsing data rows
                parts = [p.strip() for p in line.split(',')]
                if len(parts) != len(header):
                    continue
                row = dict(zip(header, parts))
                # logging.info("Debug row keys: %s", row.keys()) # Uncomment for verbose debug
                if 'S-O-T' in row:
                    try:
                        elong = float(row['S-O-T'])
                        # Gunakan kolom pertama sebagai date
                        date_col = header[0]
                        dt_str = row[date_col].strip()
                        dt = datetime.strptime(dt_str, '%Y-%b-%d %H:%M').replace(tzinfo=ZoneInfo('UTC'))
                        if elong < min_elon:
                            min_elon = elong
                            min_time = dt
                    except (ValueError, TypeError):
                        continue
        return min_time

    @staticmethod
    def azimuth_difference(az1: float, az2: float) -> float:
        """
        Menghitung selisih azimut terbatas dalam rentang (−180°, 180°].
        Nilai positif berarti objek kedua berada di timur objek pertama.
        """
        diff = az2 - az1
        while diff <= -180.0:
            diff += 360.0
        while diff > 180.0:
            diff -= 360.0
        return diff

    def compute_hilal_parameters(self,
                                 date_estimate: datetime,
                                 newmoon_search_window: timedelta = timedelta(days=2),
                                 sampling_step: timedelta = timedelta(minutes=10)) -> Optional[Dict[str, object]]:
        """
        Workflow tinggi untuk menghitung parameter hisab hilal pada hari
        konjungsi (geosentris) di lokasi pengamat ini.

        Langkah‑langkahnya sebagai berikut:

        1. Tentukan jendela pencarian sekitar estimasi tanggal hisab
           (biasanya ±1 hari). Fungsi akan mencari waktu konjungsi
           geosentris (minimum elongasi) dalam interval ini dengan
           memanggil ``find_geocentric_conjunction``.
        2. Setelah menemukan waktu konjungsi UTC (`t_conj`), tentukan
           tanggal lokal (`t_conj_local.date()`).
        3. Ambil waktu matahari terbenam pada tanggal lokal menggunakan
           ``get_sunrise_sunset`` dan ambil `sunset_utc`.
        4. Ambil ephemeris Matahari dan Bulan pada waktu `sunset_utc`.
           Data yang diambil meliputi azimut (`AZ`), ketinggian
           (`EL`), elongasi (`S-O-T`), sudut fase (`phi`), diameter
           sudut (`Ang-diam`), jarak (`delta`), dan fraksi iluminasi
           (`Illu%`).
        5. Hitung selisih azimut Bulan–Matahari.
        6. Kembalikan struktur dict berisi semua parameter.

        Parameters
        ----------
        date_estimate : datetime
            Estimasi tanggal ijtima' (misal dari kalender perhitungan).
            Gunakan zona waktu lokal pengamat untuk parameter ini.
        newmoon_search_window : timedelta, optional
            Rentang waktu di sekitar date_estimate yang akan di-scan
            untuk mencari ijtima' geosentris. Default ±2 hari (4
            total), yaitu pencarian dimulai ``date_estimate - 1 hari`` dan
            berakhir ``date_estimate + 1 hari``.
        sampling_step : timedelta, optional
            Langkah sampling untuk pencarian minimal elongasi.

        Returns
        -------
        dict or None
            Kamus berisi parameter hisab hilal. Jika terjadi
            kegagalan (misal API tidak mengembalikan data),
            mengembalikan ``None``.
        """
        # Ubah date_estimate ke UTC
        estimate_local = date_estimate.replace(tzinfo=self.timezone)
        start_utc = (estimate_local - newmoon_search_window / 2).astimezone(ZoneInfo('UTC'))
        end_utc = (estimate_local + newmoon_search_window / 2).astimezone(ZoneInfo('UTC'))
        logger.info("Mencari ijtima' geosentris antara %s dan %s", start_utc, end_utc)
        t_conj = self.find_geocentric_conjunction(start_utc, end_utc, step=sampling_step)
        if t_conj is None:
            logger.warning("Tidak menemukan waktu konjungsi dalam interval pencarian.")
            return None
        logger.info("Waktu ijtima' geosentris ditemukan: %s", t_conj)
        # Konversi ke waktu lokal untuk mengetahui tanggal kalender
        t_conj_local = t_conj.astimezone(self.timezone)
        date_local = datetime(t_conj_local.year, t_conj_local.month, t_conj_local.day, tzinfo=self.timezone)
        # Ambil waktu terbit/terbenam matahari
        sunrise_utc, sunset_utc = self.get_sunrise_sunset(date_local)
        if sunset_utc is None:
            logger.warning("Waktu matahari terbenam tidak tersedia.")
            return None
        logger.info("Matahari terbenam (UTC): %s", sunset_utc)
        # Ambil ephemeris Sun dan Moon pada saat terbenam
        times = [sunset_utc]
        # Gunakan quantities: 4 (AZ,EL) + 10 (illum), 13 (Ang-diam), 20 (delta), 23 (elong), 43 (phi)
        quantities = '4,10,13,20,23,43'
        sun_rows = self.query_ephemeris('10', times, quantities)
        moon_rows = self.query_ephemeris('301', times, quantities)
        if not sun_rows or not moon_rows:
            logger.warning("Data ephemeris tidak lengkap.")
            return None
        sun_row = sun_rows[0]
        moon_row = moon_rows[0]
        
        def get_col(row, *possible_names):
            # Helper to find value from row if key contains one of possible_names
            for key in row:
                for name in possible_names:
                    if name in key:
                         try:
                             return float(row[key])
                         except ValueError:
                             pass
            return math.nan

        try:
            sun_az = get_col(sun_row, 'Azi', 'AZ')
            sun_el = get_col(sun_row, 'Elev', 'EL')
            moon_az = get_col(moon_row, 'Azi', 'AZ')
            moon_el = get_col(moon_row, 'Elev', 'EL')
            
            elong = get_col(moon_row, 'S-O-T', 'Elong')
            phase = get_col(moon_row, 'phi', 'Phase')
            ang_diam = get_col(moon_row, 'Ang-diam', 'Diam')
            semidiameter = ang_diam / 2.0
            delta = get_col(moon_row, 'delta', 'Dist')
            illum = get_col(moon_row, 'Illu')
        except (KeyError, ValueError):
            logger.error("Gagal mengurai ephemeris.")
            return None
        # Selisih azimut Bulan – Matahari
        diff_az = self.azimuth_difference(sun_az, moon_az)
        # Jarak zenith
        zenith_distance_moon = 90.0 - moon_el
        return {
            'conjunction_utc': t_conj,
            'sunset_utc': sunset_utc,
            'sun_azimuth': sun_az,
            'sun_altitude': sun_el,
            'moon_azimuth': moon_az,
            'moon_altitude': moon_el,
            'azimuth_difference': diff_az,
            'elongation': elong,
            'phase_angle': phase,
            'illuminated_fraction': illum,
            'moon_semidiameter': semidiameter,
            'moon_range_au': delta,
            'moon_zenith_distance': zenith_distance_moon,
        }


def main():
    """
    Contoh penggunaan modul ini. Fungsi ``main`` akan membaca
    lokasi Semarang, Indonesia (110.42° BT, −7.00° LS, elevasi 0 km),
    menentukan estimasi tanggal ijtima' dari kalender (misalnya
    10 Mei 2026) dan menghasilkan parameter hisab hilal sesuai
    langkah‐langkah yang dijelaskan di ``compute_hilal_parameters``.

    Jalankan `python program_hisab_jpl_horizons.py` untuk melihat
    contoh output. Perhatikan bahwa jika dijalankan di lingkungan
    tanpa akses internet, pemanggilan API akan gagal.
    """
    # Lokasi Semarang: longitude timur positif, latitude utara positif
    location = {
        'lon': 110.42,       # derajat timur
        'lat': -6.97,        # derajat selatan → negatif
        'elevation': 0.0,    # km
    }
    client = HorizonsClient(location=location, timezone='Asia/Jakarta')
    # Estimasi tanggal ijtima' – contoh: 10 Mei 2026 (zona waktu lokal)
    estimate = datetime(2026, 5, 10, 12, 0)  # jam tengah hari lokal
    results = client.compute_hilal_parameters(estimate)
    if results:
        print("Waktu ijtima' (UTC):", results['conjunction_utc'])
        print("Waktu matahari terbenam (UTC):", results['sunset_utc'])
        print("Azimut Matahari (deg):", results['sun_azimuth'])
        print("Ketinggian Matahari (deg):", results['sun_altitude'])
        print("Azimut Bulan (deg):", results['moon_azimuth'])
        print("Ketinggian Bulan (deg):", results['moon_altitude'])
        print("Selisih azimut (Bulan–Matahari) (deg):", results['azimuth_difference'])
        print("Elongasi (deg):", results['elongation'])
        print("Sudut fase (deg):", results['phase_angle'])
        print("Fraksi iluminasi (%):", results['illuminated_fraction'])
        print("Semidiameter Bulan (arcsec):", results['moon_semidiameter'])
        print("Jarak Bulan (AU):", results['moon_range_au'])
        print("Jarak zenith Bulan (deg):", results['moon_zenith_distance'])
    else:
        print("Gagal mengambil data hilal.")


if __name__ == '__main__':
    main()