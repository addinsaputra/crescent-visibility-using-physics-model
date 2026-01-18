
import math
from visual_limit_calc2 import visual_limit
from kastner import hitung_luminansi_kastner
from schaefer_limiting_magnitude import (
    compute_limiting_magnitude, 
    SchaeferInputs, 
    eye_pupil_diameter
)

def nL_to_mag_arcsec2(nL):
    """
    Inverse of B = 34.08 * exp(20.7233 - 0.92104 * B_V)
    B/34.08 = exp(...)
    ln(B/34.08) = 20.7233 - 0.92104 * B_V
    0.92104 * B_V = 20.7233 - ln(B/34.08)
    B_V = (20.7233 - math.log(nL / 34.08)) / 0.92104
    """
    if nL <= 0:
        return 22.0 # Default dark
        
    val = nL / 34.08
    if val <= 0:
        return 22.0
    
    return (20.7233 - math.log(val)) / 0.92104

def subtend_altaz(alt1, az1, alt2, az2):
    """
    Calculate angular separation (degrees) between two points in Alt/Az coordinates.
    All inputs in degrees.
    """
    rd = math.pi / 180.0
    alt1_r = alt1 * rd
    alt2_r = alt2 * rd
    az_diff_r = (az1 - az2) * rd
    
    cos_theta = math.sin(alt1_r) * math.sin(alt2_r) + \
                math.cos(alt1_r) * math.cos(alt2_r) * math.cos(az_diff_r)
    
    # Clamp for numerical stability
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.degrees(math.acos(cos_theta))

def main():
    print("=== Perhitungan Delta m (Kastner + Schaefer) ===")
    # ---------------------------------------------------------
    # 1. Parameter Input saat matahari terbenam (Sunset)
    # ---------------------------------------------------------
    # data input schaefer
    month                = 1 # Bulan
    year                 = 2024 # Tahun
    phase_angle_schaefer = 0.0 # input tidak usah dirubah
    alt_moon             = -90.0 # input tidak usah dirubah
    az_moon_objk         = 180.0 # azimuth bulan terhadap objek
    alt_sun              = -0.27 # Altitude Matahari topocentric
    az_sun_objk          = 5.036944444 # azimuth matahari terhadap objek
    obs_humidity         = 80.66 # Kelembaban
    obs_temp             = 27.4 # Temperatur
    obs_lat              = -6.99 # Latitude
    obs_elev             = 89.0 # Elevasi
    snellen_ratio        = 1.0 # Snellen Ratio
    alt_objk             = 4.122777778 # Altitude objek topocentric

    # ---------------------------------------------------------
    # 1b. Parameter Teleskop (Faktor Koreksi)
    # ---------------------------------------------------------
    use_telescope        = True
    telescope_aperture   = 100.0  # Diameter (mm)
    telescope_focal_len  = 1000.0 # Focal Length (mm) - optional, for calc mag if needed
    # Magnifikasi: Jika user punya eyepiece tetep, masukan langsung.
    # Disini kita asumsi user input Magnifikasi
    telescope_mag        = 50.0   # Magnification (m)   
    tel_transmission     = 0.96   # e.g. 0.9^4 or similar (t)   
    tel_obstruction      = 0.0    # mm (e.g. secondary mirror) D_s/D
    tel_surfaces         = 6      # Number of optical surfaces (n)
    
    # Parameter Pengamat Tambahan
    observer_age         = 22     # Tahun (A)
    star_color_index     = 0.9    # B-V (Moon/Sun is approx 0.6-1.0)
    seeing_arcsec        = 3.0    # Seeing condition (theta_s)

    # data input kastner
    phase_angle_kastner  = 174.8791667 # fase bulan
    elongation           = 6.971944444 # Elongasi bulan terhadap objek 
    moon_semidiameter    = 0.262 # 'r' dalam derajat (approx 15-16 arcmin ~ 0.25-0.26 deg)
    zenith_dist          = 90.0 - alt_objk  # Zenith distance

    # ---------------------------------------------------------
    # 2. Hitung Sky Brightness (Schaefer)
    # ---------------------------------------------------------
    # visual_limit_calc2 return dict keys: 'sky_brightness', etc.
    res_schaefer = visual_limit(
        month=month,
        year=year,
        phase_angle=phase_angle_schaefer,
        altmoon=alt_moon,   # Bulan ada di posisi ini
        azimoon=az_moon_objk,# Relative az=0
        altsun=alt_sun,
        azisun=az_sun_objk,
        humidity=obs_humidity,
        temperature=obs_temp,
        latitude=obs_lat,
        altitude=obs_elev,
        snellen_ratio=snellen_ratio,
        altstar=alt_objk,   # Objeknya adalah bulan (Hilal)
        goal_magnitude=99   # Hitung nilai saat ini
    )
    
    sky_brightness_nL = res_schaefer["sky_brightness"]
    k_v = res_schaefer["K"][2] # Extinction coeff di V-band
    
    print(f"\n[Schaefer] Sky Brightness: {sky_brightness_nL:.4f} nL")
    print(f"[Schaefer] Extinction Coeff (k): {k_v:.4f}")

    # ---------------------------------------------------------
    # 2b. Hitung Faktor Koreksi Teleskop (Limiting Magnitude)
    # ---------------------------------------------------------
    # Kita hitung limiting magnitude untuk:
    # A. Mata Telanjang (Naked Eye)
    # B. Dengan Teleskop
    # Selisihnya adalah gain/correction factor.
    
    telescope_gain_mag = 0.0
    
    if use_telescope and sky_brightness_nL > 0:
        # Convert B(nL) to Mag/arcsec^2
        b_mag_arcsec = nL_to_mag_arcsec2(sky_brightness_nL)
        
        # 1. Naked Eye Baseline
        # Aperture = Eye Pupil, Magnification = 1
        d_eye = eye_pupil_diameter(observer_age)
        
        inputs_eye = SchaeferInputs(
            sky_brightness_mmag=b_mag_arcsec,
            aperture=d_eye,
            central_obstruction=0.0,
            transmission_per_surface=1.0, # Naked eye transmissivity considered in standard I? 
                                          # Schaefer paper implies I is eye sensitivity entering pupil. 
                                          # So transmission loss is 0 for eye entry? 
                                          # Or transmission=1.
            n_surfaces=0,
            magnification=1.0,
            seeing=seeing_arcsec,         # Seeing doesn't affect naked eye much unless extended?
            extinction_coefficient=k_v,
            zenith_distance=zenith_dist,
            age=observer_age,
            colour_index=star_color_index
        )
        
        lim_mag_eye = compute_limiting_magnitude(inputs_eye)
        
        # 2. Telescope
        inputs_scope = SchaeferInputs(
            sky_brightness_mmag=b_mag_arcsec,
            aperture=telescope_aperture,
            central_obstruction=tel_obstruction,
            transmission_per_surface=math.pow(tel_transmission, 1.0/max(1, tel_surfaces)) if tel_surfaces > 0 else 1.0, 
            # Note: inputs definition says 'transmission_per_surface'. 
            # If user supplies total transmission 'tel_transmission', we reverse it?
            # Or just assume user sets per surface?
            # Let's assume input 'tel_transmission' above is TOTAL transmission estimate?
            # Wait, best to use standard: t roughly 0.95 per surface?
            # Let's overwrite for safety to match Schaefer inputs expectations:
            # SchaeferInputs expects 'transmission_per_surface'.
            n_surfaces=tel_surfaces,
            magnification=telescope_mag,
            seeing=seeing_arcsec,
            extinction_coefficient=k_v,
            zenith_distance=zenith_dist,
            age=observer_age,
            colour_index=star_color_index
        )
        # Fix transmission inputs:
        # If we use standard 0.95 per surface:
        inputs_scope.transmission_per_surface = 0.95 
        
        lim_mag_scope = compute_limiting_magnitude(inputs_scope)
        
        telescope_gain_mag = lim_mag_scope - lim_mag_eye
        '''
        print(f"\n[Telescope Correction]")
        print(f"Limiting Mag (Naked Eye): {lim_mag_eye:.2f}")
        print(f"Limiting Mag (Telescope): {lim_mag_scope:.2f}")
        print(f"Gain (Correction Factor): {telescope_gain_mag:.2f} mag")
        '''
    # ---------------------------------------------------------
    # 3. Hitung Luminansi Hilal (Kastner)
    # ---------------------------------------------------------
    # Butuh: alpha, elongation, r, z
    
    # Hitung Elongasi (Jarak sudut Matahari - Bulan)
    
    # Gunakan k dari Schaefer untuk konsistensi? Atau default Kastner 0.2?
    # Kastner paper mungkin menyarankan k ~ 0.2, tapi Schaefer menghitung real-time.
    # Kita coba gunakan k Schaefer agar lebih akurat dengan kondisi atmosfer saat itu.
    k_used = k_v 
    
    lum_hilal_nL = hitung_luminansi_kastner(
        alpha=phase_angle_kastner,
        elongation=elongation,
        r=moon_semidiameter,
        z=zenith_dist,
        k=k_used
    )
    
    print(f"[Kastner]  Luminansi Hilal Atmosfer: {lum_hilal_nL:.4f} nL")
    
    # ---------------------------------------------------------
    # 4. Hitung Delta m
    # ---------------------------------------------------------
    # Delta_m = 2.5 * log(R)
    # R = Luminansi Hilal / Sky Brightness
    
    if sky_brightness_nL <= 0:
        print("\nError: Sky Brightness <= 0 (Terlalu gelap atau error perhitungan)")
        return

    R = lum_hilal_nL / sky_brightness_nL
    
    # Guard against log(0)
    if R <= 0:
        delta_m = -999.0 # Tidak terlihat
    else:
        delta_m = 2.5 * math.log10(R)
        
    # Delta m ini adalah naked eye (berdasarkan formula contrast ratio)
    # Jika menggunakan teleskop, visibilitas meningkat sebanding dengan gain magnitude.
    # Jadi kita tambahkan gain ke Delta m?
    # Delta m (Kastner) adalah ukuran visual. Semakin besar semakin terlihat.
    # Jika teleskop membuat kita bisa melihat +2 mag lebih redup, maka Delta m efektif naik +2?
    
    delta_m_naked = 0.0
    if R <= 0:
        delta_m_naked = -99.0
    else:
        delta_m_naked = 2.5 * math.log10(R)
        
    delta_m_final = delta_m_naked + telescope_gain_mag

    print(f"\n[HASIL] Rasio Kontras (R) Naked: {R:.6f}")
    print(f"[HASIL] Delta m (Naked Eye)    : {delta_m_naked:.4f}")
    if use_telescope:
        print(f"[HASIL] Delta m (Telescope)    : {delta_m_final:.4f}")
    else:
        delta_m_final = delta_m_naked

    # Evaluasi berdasarkan Delta m Final
    delta_m = delta_m_final
    
    if delta_m < 0:
        print("=> Kemungkinan SULIT TERLIHAT (Delta m < 0)")
    else:
        print("=> Kemungkinan TERLIHAT (Delta m > 0)")

if __name__ == "__main__":
    main()
