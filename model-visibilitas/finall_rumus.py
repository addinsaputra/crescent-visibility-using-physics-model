"""
Script untuk menghitung visibilitas hilal dengan kombinasi model Kastner
dan koreksi teleskopik Schaefer.  Program ini mengintegrasikan modul
``kastner``, ``visual_limit_calc2`` dan ``schaefer_limiting_magnitude``
yang telah diperbaiki.  Output yang dicetak mencakup luminansi hilal dan
kecerahan langit sebelum dan sesudah koreksi teleskopik, koefisien
ekstingsi, rasio kontras dan nilai Δm untuk mata telanjang dan teleskop.

Integrasi NASA POWER:
Program ini juga dapat mengambil data kelembapan relatif dan suhu udara
secara otomatis dari NASA POWER API dengan memanggil fungsi
``get_rh_t_at_time`` sebelum perhitungan visibilitas.  Data diinterpolasi
ke waktu sunset dalam UTC.

Algoritma:

1. Gunakan ``visual_limit_calc2.visual_limit`` untuk menghitung
   kecerahan langit senja (nL) dan koefisien ekstingsi ``k_V``.  Fungsi
   ini juga menghitung magnitudo limit bintang tetapi kita hanya
   membutuhkan nilai ``sky_brightness`` dan elemen ketiga dari array
   ``K`` (V‑band).

2. Gunakan ``kastner.hitung_luminansi_kastner`` dengan parameter fase,
   elongasi, semidiameter bulan, jarak zenith dan koefisien ekstingsi
   ``k_V`` untuk mendapatkan luminansi hilal di dalam atmosfer (nL).

3. Hitung rasio kontras naked eye: ``R_naked = L_hilal / B_sky`` dan
   nilai ``Δm_naked = 2.5 * log10(R_naked)``.  Nilai ∆m > 0 menunjukkan
   kontras di atas satu (hilal lebih terang daripada langit), sedangkan
   nilai negatif menunjukkan sebaliknya.

4. Bila teleskop digunakan, hitung faktor koreksi total ``F``
   menggunakan fungsi ``extended_surface_correction_factor`` dari modul
   Schaefer.  Faktor ini meliputi koreksi binokular, transmisi,
   pupil, area, Stiles–Crawford, magnifikasi (1/M²) dan warna.
   Kalikan ``L_hilal`` dan ``B_sky`` dengan ``F`` untuk mendapatkan
   luminansi hilal dan kecerahan langit yang terlihat melalui teleskop.
   Karena faktor yang sama diterapkan ke keduanya, rasio kontras tetap
   sama, tetapi besarnya pencahayaan yang diterima retina meningkat atau
   menurun sesuai dengan konfigurasi teleskop.

5. Cetak hasil sesuai dengan permintaan.

Catatan: Untuk memanfaatkan magnifikasi teleskop secara penuh, model
ambang visibilitas harus mempertimbangkan ukuran sudut hilal (misalnya
ketebalan atau panjang busur).  Program ini tidak mengubah parameter
geometri hilal karena ``kastner.hitung_luminansi_kastner`` memodelkan
kecerahan fisik, bukan ambang deteksi.  Penyesuaian terhadap ambang
deteksi (misalnya dengan mengalikan ukuran sudut dengan faktor M)
harus dilakukan pada tahap evaluasi visibilitas, bukan pada perhitungan
luminansi.
"""

import math
from datetime import datetime, timezone
from typing import Optional

from kastner import hitung_luminansi_kastner
from visual_limit_calc2 import visual_limit
from limit_teleskop import extended_surface_correction_factor
from meteo_power_api import ObservingLocation, get_rh_t_at_time, PowerAPIError


def compute_hilal_visibility(
    month: int,
    year: int,
    phase_angle_schaefer: float,
    altmoon: float,
    azimoon: float,
    altsun: float,
    azisun: float,
    humidity: float,
    temperature: float,
    latitude: float,
    altitude: float,
    snellen_ratio: float,
    altstar: float,
    phase_angle_kastner: float,
    elongation: float,
    moon_semidiameter: float,
    zenith_distance: float,
    use_telescope: bool = True,
    aperture: float = 100.0,
    magnification: float = 50.0,
    transmission: float = 0.95,
    n_surfaces: int = 6,
    central_obstruction: float = 0.0,
    observer_age: float = 30.0,
    colour_index: float = 0.9,
) -> None:
    """Hitung dan cetak parameter visibilitas hilal.

    Parameter yang diterima sama seperti pada program contoh yang ada
    sebelumnya.  Jika teleskop digunakan, maka parameter ``aperture``
    (mm), ``magnification``, ``transmission`` (total atau per surface),
    ``n_surfaces`` (jumlah permukaan optik) dan ``central_obstruction``
    (diameter obstruction) akan digunakan untuk menghitung faktor
    koreksi teleskopik.
    """
    # ------------------------------------------------------------------
    # 1. Hitung kecerahan langit dan koefisien ekstingsi menggunakan
    #    model Schaefer/Reijs.  Keluaran fungsi adalah dictionary.
    #    Jika goal_magnitude=99 (atau None) maka program menghitung
    #    magnitudo limit untuk ketinggian objek yang diberikan.
    result = visual_limit(
        month=month,
        year=year,
        phase_angle=phase_angle_schaefer,
        altmoon=altmoon,
        azimoon=azimoon,
        altsun=altsun,
        azisun=azisun,
        humidity=humidity,
        temperature=temperature,
        latitude=latitude,
        altitude=altitude,
        snellen_ratio=snellen_ratio,
        altstar=altstar,
        goal_magnitude=99,
    )
    # Sky brightness in nanoLambert (nL) – note: 1 nL = 1 mμL
    sky_brightness_nl = float(result.get("sky_brightness", 0.0))
    # Extinction coefficient per airmass; take V-band (index 2) from K.
    # Fallback to DM if K is missing for compatibility.
    K_list = result.get("K", [])
    DM_list = result.get("DM", [])
    k_v: float = 0.0
    if len(K_list) > 2:
        k_v = float(K_list[2])
    elif len(DM_list) > 2:
        k_v = float(DM_list[2])

    # ------------------------------------------------------------------
    # 2. Hitung luminansi hilal menggunakan model Kastner.  Parameter
    #    k diambil dari koefisien ekstingsi yang dihitung di atas.
    lum_hilal_nl = hitung_luminansi_kastner(
        alpha=phase_angle_kastner,
        elongation=elongation,
        r=moon_semidiameter,
        z=zenith_distance,
        k=k_v,
    )

    # ------------------------------------------------------------------
    # 3. Hitung rasio kontras dan Δm untuk mata telanjang.
    if sky_brightness_nl <= 0.0:
        raise ValueError(
            "Kecerahan langit harus positif untuk menghitung rasio kontras."
        )
    ratio_naked = lum_hilal_nl / sky_brightness_nl
    delta_m_naked: Optional[float]
    if ratio_naked > 0.0:
        delta_m_naked = 2.5 * math.log10(ratio_naked)
    else:
        delta_m_naked = None

    # ------------------------------------------------------------------
    # 4. Bila teleskop digunakan, hitung faktor koreksi.
    if use_telescope:
        # Perkirakan transmisi per permukaan.
        if n_surfaces > 0 and transmission > 0:
            t_per_surface = transmission ** (1.0 / n_surfaces)
        else:
            t_per_surface = transmission
            
        # Hitung faktor total untuk EXTENDED source (Langit)
        # Faktor ini mencakup dimming akibat magnifikasi (1/M^2)
        F_sky = extended_surface_correction_factor(
            B=sky_brightness_nl,
            aperture=aperture,
            magnification=magnification,
            central_obstruction=central_obstruction,
            transmission_per_surface=t_per_surface,
            n_surfaces=n_surfaces,
            age=observer_age,
            colour_index=colour_index,
        )
        
        # Hitung faktor total untuk HILAL (Signal Source)
        # Hilal dianggap sebagai "signal" yang menderita loss transmisi dsb,
        # TAPI magnifikasi mendistribusikan cahaya ke retina yang lebih luas 
        # (atau untuk point source, tidak meredupkan surface brightness objeknya).
        # Untuk mendapatkan "Contrast Gain", kita asumsikan Hilal tidak mengalami 
        # dimming 1/M^2 dari perspektif kontras terhadap background.
        # Jadi kita kalikan F_sky dengan M^2 untuk membatalkan efek 1/M^2.
        F_hilal = F_sky * (magnification ** 2)

        sky_brightness_tel = sky_brightness_nl * F_sky
        lum_hilal_tel = lum_hilal_nl * F_hilal
        
        # Rasio teleskop sekarang berbeda karena F_hilal != F_sky
        ratio_tel = lum_hilal_tel / sky_brightness_tel if sky_brightness_tel > 0 else float("nan")
        
        delta_m_tel: Optional[float]
        if ratio_tel > 0.0:
            delta_m_tel = 2.5 * math.log10(ratio_tel)
        else:
            delta_m_tel = None
    else:
        F_sky = 1.0
        F_hilal = 1.0
        sky_brightness_tel = sky_brightness_nl
        lum_hilal_tel = lum_hilal_nl
        ratio_tel = ratio_naked
        delta_m_tel = delta_m_naked

    # ------------------------------------------------------------------
    # 5. Cetak hasil
    print("=== Hasil Perhitungan Visibilitas Hilal ===")
    print(f"Luminansi Hilal (nL) Naked Eye       : {lum_hilal_nl}")
    print(f"Sky Brightness (nL) Naked Eye        : {sky_brightness_nl}")
    if use_telescope:
        print(f"Luminansi Hilal (nL) Teleskop    : {lum_hilal_tel}")
        print(f"Sky Brightness (nL) Teleskop     : {sky_brightness_tel}")
    print(f"Koefisien Extingsi (k_V)             : {k_v}")
    print(f"Rasio Kontras R Naked Eye            : {ratio_naked}")
    if use_telescope:
        print(f"Rasio Kontras R Teleskop         : {ratio_tel}")
    # Tampilkan Δm jika terdefinisi; gunakan simbol '-' untuk tidak terdefinisi
    def fmt_delta(dm: Optional[float]) -> str:
        return f"{dm}" if dm is not None else "-"

    print(f"Delta m (Naked Eye)                : {fmt_delta(delta_m_naked)}")
    if use_telescope:
        print(f"Delta m (Teleskop)                 : {fmt_delta(delta_m_tel)}")


def compute_hilal_visibility_with_power(
    location_name: str,
    latitude: float,
    longitude: float,
    altitude: float,
    timezone_str: str,
    sunset_utc: datetime,
    month: int,
    year: int,
    phase_angle_schaefer: float,
    altmoon: float,
    azimoon: float,
    altsun: float,
    azisun: float,
    snellen_ratio: float,
    altstar: float,
    phase_angle_kastner: float,
    elongation: float,
    moon_semidiameter: float,
    zenith_distance: float,
    use_telescope: bool = True,
    aperture: float = 100.0,
    magnification: float = 50.0,
    transmission: float = 0.95,
    n_surfaces: int = 6,
    central_obstruction: float = 0.0,
    observer_age: float = 30.0,
    colour_index: float = 0.9,
) -> None:
    """Hitung visibilitas hilal dengan data RH dan suhu dari NASA POWER.

    Fungsi ini mengambil data kelembapan relatif dan suhu udara dari
    NASA POWER API, menginterpolasinya ke waktu sunset UTC, lalu
    memanggil ``compute_hilal_visibility`` dengan nilai tersebut.

    Parameters
    ----------
    location_name : str
        Nama lokasi pengamatan.
    latitude : float
        Lintang dalam derajat desimal.
    longitude : float
        Bujur dalam derajat desimal.
    altitude : float
        Ketinggian dalam meter di atas permukaan laut.
    timezone_str : str
        Zona waktu IANA (misal "Asia/Jakarta").
    sunset_utc : datetime
        Waktu sunset dalam UTC (harus timezone-aware).
    [parameter lainnya sama dengan compute_hilal_visibility]
    """
    loc = ObservingLocation(
        name=location_name,
        latitude=latitude,
        longitude=longitude,
        altitude=altitude,
        timezone=timezone_str,
    )

    print(f"Mengambil data dari NASA POWER untuk {location_name}...")
    print(f"Waktu sunset UTC: {sunset_utc.isoformat()}")

    try:
        rh, temp = get_rh_t_at_time(loc, sunset_utc)
        print(f"NASA POWER - Kelembapan Relatif: {rh}%")
        print(f"NASA POWER - Suhu Udara: {temp}°C")
        print("-" * 50)
    except PowerAPIError as e:
        print(f"Error mengambil data NASA POWER: {e}")
        raise

    compute_hilal_visibility(
        month=month,
        year=year,
        phase_angle_schaefer=phase_angle_schaefer,
        altmoon=altmoon,
        azimoon=azimoon,
        altsun=altsun,
        azisun=azisun,
        humidity=rh,
        temperature=temp,
        latitude=latitude,
        altitude=altitude,
        snellen_ratio=snellen_ratio,
        altstar=altstar,
        phase_angle_kastner=phase_angle_kastner,
        elongation=elongation,
        moon_semidiameter=moon_semidiameter,
        zenith_distance=zenith_distance,
        use_telescope=use_telescope,
        aperture=aperture,
        magnification=magnification,
        transmission=transmission,
        n_surfaces=n_surfaces,
        central_obstruction=central_obstruction,
        observer_age=observer_age,
        colour_index=colour_index,
    )


if __name__ == "__main__":
    # Contoh penggunaan MANUAL tanpa NASA POWER API
    # Data RH dan Temperature dimasukkan secara manual
    
    compute_hilal_visibility(
        month=7,
        year=2022,
        phase_angle_schaefer=0.0,
        altmoon=-90.0,
        azimoon=180.0,
        altsun=-0.3,
        azisun=3.037,
        humidity=89.0,  # Manual input
        temperature=26.0,   # Manual input
        latitude=-6.20889,
        altitude=89.0,
        snellen_ratio=1.0,
        altstar=6.638611111,
        phase_angle_kastner=171.5919444,
        elongation=8.385833333,
        moon_semidiameter=0.2625,
        zenith_distance=90.0 - 6.638611111,
        use_telescope=True,
        aperture=100.0,
        magnification=50.0,
        transmission=0.96,
        n_surfaces=6,
        central_obstruction=0.0,
        observer_age=22.0,
        colour_index=0.9,
    )
