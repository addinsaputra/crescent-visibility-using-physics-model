
import math

"""
This module combines the Kastner visibility model for the lunar crescent
with a physically motivated implementation of the telescope correction
factors described in Bradley E. Schaefer’s 1990 paper “Telescopic
Limiting Magnitudes.”  The original version of this script estimated
how much a telescope improves visibility by calling a limiting magnitude
solver for point sources and simply adding the magnitude gain onto
Kastner’s contrast ratio.  That approach is not appropriate for the
lunar crescent because Schaefer’s limiting‑magnitude formulae apply
to point sources (stars)【487407022821281†L82-L100】 and explicitly note that
their extension to extended sources is a difficult problem【487407022821281†L45-L47】.  In
particular, the gain in limiting magnitude for stars arises from
concentrating the light of a point source, whereas a crescent moon
and the twilight sky are extended sources whose surface brightnesses
are equally affected by the telescope optics【487407022821281†L186-L193】.  Simply
adding a star‑magnitude gain to an extended‑source contrast will
overestimate the improvement and, in some cases, change the contrast
in the wrong direction.

To remedy this, the revised code below implements the individual
correction factors from Schaefer’s eqs. (14)–(15) explicitly and
applies them consistently to both the hilal (extended source) and
the twilight sky.  For extended sources the magnification factor
must be applied to the object as well as the background【487407022821281†L186-L193】.
The overall effect is that the surface brightnesses of the crescent
and sky are both scaled by the same factor, so their ratio does not
change; the principal benefit of a telescope is the increase in the
apparent angular size of the crescent, which reduces the contrast
threshold.  Implementing a proper size‑dependent threshold lies
beyond the scope of this example, but the code prepares the
brightnesses for such a model.
"""

from visual_limit_calc2 import visual_limit
from kastner import hitung_luminansi_kastner

try:
    # We import eye_pupil_diameter from the Schaefer implementation
    # if available; otherwise we provide a local definition below.
    from schaefer_limiting_magnitude import eye_pupil_diameter  # type: ignore
except Exception:
    def eye_pupil_diameter(age: float) -> float:
        """Return the eye pupil diameter in mm as a function of age.

        Schaefer (1990) combines pupil‑diameter data from Kumnick (1954)
        and Kadlecova et al. (1958) into a convenient analytic form:
        D_e = 7 mm exp(−0.5 (A/100)^2)【487407022821281†L170-L176】.
        """
        return 7.0 * math.exp(-0.5 * ((age) / 100.0) ** 2)

def nL_to_mag_arcsec2(nL: float) -> float:
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
    
    # ---------------------------------------------------------------------
    # 2b. Hitung Faktor Koreksi Teleskop untuk sumber terentang (extended)
    # ---------------------------------------------------------------------
    # Untuk hilal dan langit senja (keduanya sumber terentang), persamaan (15)
    # Schaefer menyatakan bahwa brightness yang dilihat melalui teleskop
    # diperoleh dengan mengalikan brightness asli dengan produk faktor:
    # Fb, Ft, Fp, Fa, Fsc, Fm, Fc【487407022821281†L328-L333】.  Fm=1/M^2
    # hanya disebut diterapkan pada background dalam paper karena target
    # yang dibahas adalah point source, tetapi untuk extended source seperti
    # hilal dan langit, faktor ini berlaku pada keduanya【487407022821281†L186-L193】.
    telescope_gain_mag = 0.0  # diinisialisasi tetapi tidak lagi digunakan
    
    def compute_telescope_correction_extended(
        B_nL: float,
        D_mm: float,
        M: float,
        Ds_mm: float,
        transmission_per_surface: float,
        n_surfaces: int,
        age: float,
        colour_index: float
    ) -> float:
        """Compute the combined telescope correction factor for an extended source.

        Parameters
        ----------
        B_nL : float
            Observed sky brightness in nanoLambert (nL).
        D_mm : float
            Telescope aperture diameter in millimetres.
        M : float
            Magnification of the telescope.
        Ds_mm : float
            Central obstruction diameter in millimetres (0 for unobstructed).
        transmission_per_surface : float
            Transmission per optical surface (0–1).  Typical coated optics are
            around 0.95–0.97 per surface【487407022821281†L149-L164】.
        n_surfaces : int
            Total number of optical surfaces that light passes through.
        age : float
            Age of the observer in years.
        colour_index : float
            Colour index (B−V) of the source; used for colour correction.

        Returns
        -------
        float
            Multiplicative factor by which both the source and background
            brightnesses should be scaled when viewed through the telescope.

        Notes
        -----
        The factor returned is

        Fb · Ft · Fp · Fa · Fsc · Fm · Fc

        following eq. (15) of Schaefer (1990)【487407022821281†L328-L333】, with the
        magnification factor Fm applied to extended sources as discussed
        in the text.
        """
        # Guard against non‑positive brightness
        if B_nL <= 0.0:
            return 1.0
        # Binocular to monocular correction (Fb).  Knoll et al. used
        # binocular vision, whereas telescopic observing is monocular【487407022821281†L110-L118】.
        F_b = 1.41
        # Telescope transmission and central obstruction correction (Ft).
        D_ratio = Ds_mm / D_mm if D_mm > 0 else 0.0
        Ft = 1.0
        if D_mm > 0.0:
            aperture_loss = (1.0 - D_ratio ** 2)
            # Avoid division by zero if the obstruction equals the aperture
            if aperture_loss > 0.0 and transmission_per_surface > 0.0:
                Ft = 1.0 / ((transmission_per_surface ** n_surfaces) * aperture_loss)
        # Pupil/exiting pupil loss (Fp).  If the telescope’s exit pupil
        # is larger than the eye’s pupil, some light is lost【487407022821281†L167-L179】.
        D_e = eye_pupil_diameter(age)
        exit_pupil = D_mm / M if M > 0 else float('inf')
        F_p = 1.0
        if D_e < exit_pupil and D_e > 0:
            # Too much light exits, so we lose (D/(M D_e))^2
            F_p = (D_mm / (M * D_e)) ** 2
        # Aperture collecting area ratio (Fa).  A telescope collects
        # (D/D_e)^2 times more light than the eye【487407022821281†L183-L186】.
        F_a = 1.0
        if D_e > 0:
            F_a = (D_mm / D_e) ** 2
        # Magnification reduces the surface brightness of an extended source
        # by M^2【487407022821281†L186-L193】.
        F_m = 1.0 / (M ** 2) if M > 0 else 1.0
        # Stiles–Crawford effect (Fsc).  Depends on background brightness and
        # pupil utilisation【487407022821281†L229-L252】.
        # Convert brightness to log10(B) in nL units.  If B is very low,
        # use a floor to avoid math domain errors.
        log_B = math.log10(B_nL) if B_nL > 0 else -math.inf
        F_sc = 1.0
        # Only apply if the exit pupil is smaller than the eye’s pupil; otherwise
        # the telescope fills the pupil and no edge‑falloff occurs.
        if D_e > exit_pupil and D_mm > 0 and M > 0:
            if log_B > 3.17:
                # Photopic regime (day vision).
                numerator = (D_e * M / D_mm) * (1.0 - math.exp(-0.026 * (D_mm / M) ** 2))
                denominator = (1.0 - math.exp(-0.026 * (D_mm) ** 2))
                if denominator != 0.0:
                    F_sc = numerator / denominator
            else:
                # Scotopic regime (night vision).
                numerator = 1.0 - (D_mm / (12.4 * M)) ** 4
                denominator = 1.0 - (D_e / 12.4) ** 4
                if denominator != 0.0:
                    F_sc = numerator / denominator
        # Colour correction (Fc).  For night vision the reported intensities
        # (V‑band magnitudes) must be scaled to account for the different
        # spectral response of the eye【487407022821281†L286-L317】.  For bright
        # backgrounds (log B > 3.17) no colour correction is needed.
        F_c = 1.0
        if log_B < 3.17:
            # Equation (13): −2.5 log(Fc) = 1 − (B−V)/2
            exponent = - (1.0 - (colour_index) / 2.0) / 2.5
            F_c = 10.0 ** (exponent)
        # Multiply all factors.
        F_total = F_b * Ft * F_p * F_a * F_sc * F_m * F_c
        return F_total
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
    # 4. Hitung Delta m berdasarkan rasio luminansi
    # ---------------------------------------------------------
    # Rasio kontras untuk mata telanjang (naked‑eye)
    if sky_brightness_nL <= 0.0:
        print("\nError: Sky Brightness <= 0 (Terlalu gelap atau error perhitungan)")
        return

    R_naked = lum_hilal_nL / sky_brightness_nL
    if R_naked <= 0.0:
        delta_m_naked = -99.0
    else:
        delta_m_naked = 2.5 * math.log10(R_naked)

    # Bila teleskop digunakan, terapkan faktor‑faktor koreksi untuk extended source.
    if use_telescope:
        # Convert central obstruction from mm to mm ratio; input variable
        # tel_obstruction is a linear diameter (mm).  We convert to a ratio
        # relative to the aperture to feed into the correction function.
        Ds_ratio = 0.0
        if telescope_aperture > 0.0:
            Ds_ratio = tel_obstruction / telescope_aperture
        # Use user‑supplied per‑surface transmission if provided; if the user
        # supplies total transmission (tel_transmission) we treat it as per
        # surface when n_surfaces=1.  A typical coated optical surface has
        # t≈0.95【487407022821281†L149-L164】.
        transmission_per_surface = tel_transmission
        if tel_surfaces > 0 and tel_transmission > 0:
            # Compute approximate per‑surface transmission if the user
            # specified a total transmission for all surfaces.
            transmission_per_surface = tel_transmission ** (1.0 / tel_surfaces)
        # Compute the combined correction factor.
        F_total = compute_telescope_correction_extended(
            sky_brightness_nL,
            telescope_aperture,
            telescope_mag,
            tel_obstruction,
            transmission_per_surface,
            tel_surfaces,
            observer_age,
            star_color_index
        )
        # Apply the same factor to both sky and hilal brightnesses.  This
        # preserves their ratio but yields the *perceived* luminances through
        # the telescope.
        B_tel_nL = sky_brightness_nL * F_total
        L_tel_nL = lum_hilal_nL * F_total
        # For extended sources the contrast ratio is unchanged by the
        # telescope optics【487407022821281†L186-L193】.  We compute it for completeness.
        R_tel = L_tel_nL / B_tel_nL if B_tel_nL > 0.0 else 0.0
        if R_tel <= 0.0:
            delta_m_tel = -99.0
        else:
            delta_m_tel = 2.5 * math.log10(R_tel)
        print(f"\n[HASIL] Rasio Kontras (R) Naked   : {R_naked:.6f}")
        print(f"[HASIL] Delta m (Naked Eye)      : {delta_m_naked:.4f}")
        print(f"[HASIL] Rasio Kontras (R) Telescope: {R_tel:.6f}")
        print(f"[HASIL] Delta m (Telescope)      : {delta_m_tel:.4f}")
        # Interpret the result.  Delta m > 0 implies contrast above unity.
        if delta_m_tel < 0.0:
            print("=> Kemungkinan SULIT TERLIHAT dengan teleskop (Delta m < 0)")
        else:
            print("=> Kemungkinan TERLIHAT dengan teleskop (Delta m > 0)")
    else:
        # No telescope: use naked‑eye contrast only.
        print(f"\n[HASIL] Rasio Kontras (R) Naked: {R_naked:.6f}")
        print(f"[HASIL] Delta m (Naked Eye)    : {delta_m_naked:.4f}")
        if delta_m_naked < 0.0:
            print("=> Kemungkinan SULIT TERLIHAT (Delta m < 0)")
        else:
            print("=> Kemungkinan TERLIHAT (Delta m > 0)")

if __name__ == "__main__":
    main()
