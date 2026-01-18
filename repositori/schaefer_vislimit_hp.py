"""
High‑precision implementation of Bradley E. Schaefer’s visual limiting magnitude
algorithm using modern astronomical libraries.

This module exposes a function ``vislimit_hp`` that replicates the behaviour of
Schaefer’s original VISLIMIT program but replaces the low‑precision Sun and
Moon ephemerides with high‑precision coordinates computed from NASA JPL
Development Ephemeris files via the Skyfield library.  Skyfield relies on
professional ephemeris data to deliver positions of planets and their moons
with high accuracy over centuries【708606904277315†L2-L9】.  The latest short
files in the DE series incorporate the most accurate observational data from
telescopes and spacecraft【708606904277315†L124-L137】.

The remaining steps—airmass, extinction, sky brightness, phase functions,
pupil size, critical visual angle and threshold of visibility—follow the
formulas documented in the original Schaefer algorithm and its JavaScript
translation.  Only the V band (0.55 µm) is considered, as it best matches
human photopic sensitivity【665925694870024†L455-L486】.

Usage:

    from schaefer_vislimit_hp import ObservingConditions, vislimit_hp

    obs = ObservingConditions(
        latitude=–6.966667,    # degrees (e.g. Semarang, Indonesia)
        longitude=110.416667, # east positive (degrees)
        elevation=10.0,       # metres above sea level
        humidity=70.0,        # percent
        temperature=25.0,     # °C
        age=35.0,             # years
        acuity=1.0,           # normal visual acuity
        snellen_factor=6.0    # Snellen 6/6 (20/20) eyesight
    )
    # Compute limiting magnitude at 2025‑12‑23 18:00 UTC for a star at 15° altitude
    # and azimuth 250°.
    maglim = vislimit_hp(obs, (2025, 12, 23), 18.0, 15.0, 250.0)
    print(f"Limiting magnitude: {maglim:.2f}")

Note:
    • This module requires the ``skyfield`` package and its dependency ``jplephem``
      to be installed.  Skyfield can be installed with ``pip install skyfield``
      and will download the required ephemeris file on the first run.
    • The function will raise an ImportError if Skyfield is unavailable.
    • The ephemeris used defaults to the DE440 short file which covers 1849–2150
      and offers high accuracy for planets and the Moon【708606904277315†L124-L137】.

"""

from __future__ import annotations

import math
import datetime as _datetime
from dataclasses import dataclass
from typing import Tuple

try:
    from skyfield.api import load, wgs84  # type: ignore
    from skyfield.api import Star
    _SKYFIELD_AVAILABLE = True
except ImportError:
    # Skyfield is not available in this environment.  The functions will raise
    # informative errors if called without Skyfield installed.
    _SKYFIELD_AVAILABLE = False


def julian_date(year: int, month: int, day: int, hour: float) -> float:
    """Compute the Julian Date for a calendar date and fractional hour.

    This function is provided for compatibility with Schaefer’s original
    algorithm.  It is used to compute the local sidereal time for converting
    the input star altitude and azimuth to right ascension and declination.
    """
    y = year
    m = month
    d = day + hour / 24.0
    if m <= 2:
        y -= 1
        m += 12
    a = y // 100
    b = 2 - a + a // 4
    jd = math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + d + b - 1524.5
    return jd


def local_sidereal_time(jd: float, longitude_deg: float) -> float:
    """Compute the local mean sidereal time (in hours) at the given JD and longitude.

    This formula follows the Astronomical Almanac and is identical to the one
    used in the original VISLIMIT program【665925694870024†L267-L299】.  It is used
    here for converting the star’s altitude/azimuth to equatorial coordinates.
    """
    jd_int = math.floor(jd)
    jd_frac = jd - jd_int
    if jd_frac < 0.5:
        jd_mid = jd_int - 0.5
        ut = jd_frac + 0.5
    else:
        jd_mid = jd_int + 0.5
        ut = jd_frac - 0.5
    t = (jd_mid - 2451545.0) / 36525.0
    sid_g = (24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 6.2e-6 * t * t * t) / 86400.0
    sid_g -= math.floor(sid_g)
    sid_g += 1.0027379093 * ut - longitude_deg / 360.0
    sid_g = (sid_g - math.floor(sid_g)) * 24.0
    if sid_g < 0:
        sid_g += 24.0
    return sid_g


def altaz_to_radec(
    alt_deg: float, az_deg: float, lat_deg: float, lst_hours: float
) -> Tuple[float, float]:
    """Convert an altitude/azimuth pair to equatorial coordinates.

    Given an object at altitude ``alt_deg`` and azimuth ``az_deg`` (measured
    eastwards from north) observed from latitude ``lat_deg`` with local
    sidereal time ``lst_hours``, return its right ascension (in hours) and
    declination (in degrees).  This uses the same conversion as the
    original VISLIMIT program.
    """
    alt_rad = math.radians(alt_deg)
    az_rad = math.radians(az_deg)
    lat_rad = math.radians(lat_deg)
    sin_dec = math.sin(alt_rad) * math.sin(lat_rad) + math.cos(alt_rad) * math.cos(lat_rad) * math.cos(az_rad)
    dec_rad = math.asin(sin_dec)
    ha_rad = math.atan2(
        -math.sin(az_rad) * math.cos(alt_rad),
        math.cos(az_rad) * math.cos(alt_rad) * math.sin(lat_rad) - math.tan(alt_rad) * math.cos(lat_rad),
    )
    ha_hours = ha_rad * (12.0 / math.pi)
    ra_hours = (lst_hours - ha_hours) % 24.0
    dec_deg = math.degrees(dec_rad)
    return ra_hours, dec_deg


def subtend(ra1_h: float, dec1_d: float, ra2_h: float, dec2_d: float) -> float:
    """Angular separation between two equatorial positions in radians.

    The arguments are right ascensions in hours and declinations in degrees.
    Returns the great‑circle angular distance using the spherical law of
    cosines.  Clamp the cosine to the [–1, 1] interval to avoid numerical
    errors.
    """
    ra1 = ra1_h * (math.pi / 12.0)
    dec1 = math.radians(dec1_d)
    ra2 = ra2_h * (math.pi / 12.0)
    dec2 = math.radians(dec2_d)
    cos_theta = (
        math.sin(dec1) * math.sin(dec2)
        + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    )
    cos_theta = max(-1.0, min(1.0, cos_theta))
    return math.acos(cos_theta)


@dataclass
class ObservingConditions:
    """Container for observer and environmental parameters.

    Attributes:
        latitude: Geographic latitude in degrees (north positive).
        longitude: Geographic longitude in degrees (east positive).
        elevation: Height above sea level in metres.
        humidity: Relative humidity in percent (0–100).
        temperature: Air temperature in degrees Celsius.
        age: Observer’s age in years.
        acuity: Visual acuity factor (1 = normal).
        snellen_factor: Snellen n/r ratio; 6 for 6/6 (20/20) vision.
    """
    latitude: float
    longitude: float
    elevation: float
    humidity: float
    temperature: float
    age: float
    acuity: float
    snellen_factor: float


def _ensure_skyfield():
    """Raise an informative error if Skyfield is not available."""
    if not _SKYFIELD_AVAILABLE:
        raise ImportError(
            "Skyfield is required for vislimit_hp but is not installed. "
            "Please install skyfield with 'pip install skyfield' and retry."
        )


def vislimit_hp(
    obs: ObservingConditions,
    date: Tuple[int, int, int],
    ut_hour: float,
    alt_star: float,
    az_star: float,
    ephem: str = "de440s.bsp",
) -> float:
    """Compute the visual limiting magnitude using high‑precision ephemerides.

    This function follows the same logic as the original Schaefer algorithm
    but uses Skyfield to compute the topocentric positions of the Sun and Moon
    with JPL Development Ephemerides【708606904277315†L2-L9】.  It returns the
    limiting magnitude of a point source in the V band.

    Args:
        obs: ObservingConditions describing the observer and environment.
        date: Tuple of (year, month, day).
        ut_hour: Universal Time in decimal hours (e.g. 18.5 for 18:30 UT).
        alt_star: Altitude of the target star in degrees (0–90).
        az_star: Azimuth of the target star in degrees east of north (0–360).
        ephem: Name of the JPL ephemeris file to load (defaults to DE440s).

    Returns:
        Limiting visual magnitude in the V band.

    Raises:
        ImportError: If Skyfield is not installed.
    """
    _ensure_skyfield()
    alt_star = max(0.0, min(90.0, alt_star))
    year, month, day = date
    # Convert UT hour into datetime components
    frac_hours, int_hours = math.modf(ut_hour)
    hour = int(int_hours)
    minute = int(frac_hours * 60.0)
    second = (frac_hours * 60.0 - minute) * 60.0
    dt = _datetime.datetime(year, month, day, hour, minute, int(second), int((second - int(second)) * 1e6), tzinfo=_datetime.timezone.utc)
    ts = load.timescale()
    t = ts.utc(dt)
    # Load JPL ephemeris
    planets = load(ephem)
    earth = planets["earth"]
    sun = planets["sun"]
    moon = planets["moon"]
    # Observer location
    site = wgs84.latlon(obs.latitude, obs.longitude, obs.elevation)
    observer = earth + site
    sun_ap = observer.at(t).observe(sun).apparent()
    moon_ap = observer.at(t).observe(moon).apparent()
    # Right ascension, declination and altitudes
    ra_sun = sun_ap.radec()[0].hours
    dec_sun = sun_ap.radec()[1].degrees
    sun_alt, sun_az, _sun_dist = sun_ap.altaz(
        temperature_C=obs.temperature, pressure_mbar=1010.0
    )
    alt_sun = sun_alt.degrees
    ra_moon = moon_ap.radec()[0].hours
    dec_moon = moon_ap.radec()[1].degrees
    moon_alt, moon_az, moon_dist = moon_ap.altaz(
        temperature_C=obs.temperature, pressure_mbar=1010.0
    )
    alt_moon = moon_alt.degrees
    dist_moon = moon_dist.km / 6378.14  # convert km to Earth radii
    # Compute LST using Skyfield gmst
    lst = (t.gmst + obs.longitude / 15.0) % 24.0
    # Convert star alt/az to RA/Dec
    ra_star, dec_star = altaz_to_radec(alt_star, az_star, obs.latitude, lst)
    # Zenith distances
    zz_star = math.radians(90.0 - alt_star)
    zz_moon = math.radians(90.0 - alt_moon)
    zz_sun = math.radians(90.0 - alt_sun)
    # Airmass for star
    xg = 1.0 / (math.cos(zz_star) + 0.0286 * math.exp(-10.5 * math.cos(zz_star)))
    xa = 1.0 / (math.cos(zz_star) + 0.0123 * math.exp(-24.5 * math.cos(zz_star)))
    xo = 1.0 / math.sqrt(1.0 - (math.sin(zz_star) / (1.0 + 20.0 / 6378.0)) ** 2)
    wa = 0.55
    oz = 0.031
    wt = 0.031
    kr = 0.1066 * math.exp(-obs.elevation / 8200.0) * (wa / 0.55) ** (-4.0)
    humid_frac = max(0.01, obs.humidity / 100.0)
    ka = (
        0.12
        * (wa / 0.55) ** (-1.3)
        * math.exp(-obs.elevation / 1500.0)
        * (1.0 - 0.32 / math.log(humid_frac)) ** (4.0 / 3.0)
        * (1.0 + 0.33 * math.sin(ra_sun * (math.pi / 12.0)))
    )
    lat_rad = math.radians(obs.latitude)
    ko = oz * (3.0 + 0.4 * (lat_rad * math.cos(ra_sun * (math.pi / 12.0)) - math.cos(3.0 * lat_rad))) / 3.0
    kw = wt * 0.94 * humid_frac * math.exp(obs.temperature / 15.0) * math.exp(-obs.elevation / 8200.0)
    dm_tot = kr * xg + ka * xa + ko * xo + kw * xg
    # Airmass function for Sun, Moon and star
    def airmass(zr: float) -> float:
        am = 1.0 / (math.cos(zr) + 0.025 * math.exp(-11.0 * math.cos(zr)))
        return 40.0 if zr > math.pi / 2 else am
    am_star = airmass(zz_star)
    am_moon = airmass(zz_moon)
    am_sun = airmass(zz_sun)
    # Angular separations
    rho_moon = subtend(ra_star, dec_star, ra_moon, dec_moon)
    rho_sun = subtend(ra_star, dec_star, ra_sun, dec_sun)
    rho_sm = subtend(ra_moon, dec_moon, ra_sun, dec_sun)
    f_illum = 50.0 * (1.0 + math.cos(rho_sm))
    qphase = math.acos(2.0 * f_illum / 100.0 - 1.0)
    def phase_function(rho: float) -> float:
        return 10.0 ** 5.36 * (1.06 + math.cos(rho) ** 2)
    ff_moon = phase_function(rho_moon)
    ff_sun = phase_function(rho_sun)
    def moon_brightness(f_illum: float, qphase: float) -> float:
        return -1.805 - 2.135 * math.log10(f_illum / 100.0) + 0.032 * (qphase * 180.0 / math.pi)
    smoon = moon_brightness(f_illum, qphase)
    ci_moon = ff_moon * 10 ** (-0.4 * smoon) * 10 ** (-0.4 * dm_tot * am_moon)
    ci_sun = ff_sun * 10 ** (-0.4 * (-26.74)) * 10 ** (-0.4 * dm_tot * am_sun)
    snight = 10 ** ((-22.0 + 0.28 * obs.humidity / 100.0) / 2.5)
    if alt_sun <= -18.0:
        stwilight = 0.0
    elif alt_sun >= -12.0:
        stwilight = 10 ** ((-4.0 * (alt_sun + 12.0) - 8.0) / 2.5)
    else:
        stwilight = 10 ** (-(14.4 + 0.2 * alt_sun) / 2.5)
    sky_brightness = snight + stwilight + ci_moon + ci_sun
    mlim_eff = -2.5 * math.log10(sky_brightness)
    pupil = 7.0 - 0.1 * (obs.age - 30.0)
    pupil = max(2.0, min(7.0, pupil))
    critical_angle = 0.17 * (pupil / 7.0) ** 0.4
    star_factor = 1.0
    if alt_sun > -9.0:
        c1 = 2.58
        c2 = 1.20
    else:
        c1 = 1.45
        c2 = 0.93
    ext_corr = dm_tot * am_star
    m_threshold = mlim_eff + c1 + c2 * math.log10(star_factor) + ext_corr
    acuity_corr = -2.5 * math.log10(obs.acuity) - 2.5 * math.log10(obs.snellen_factor / 6.0)

    # --- Print Outputs for Terminal ---
    k_zenith = kr + ka + ko + kw
    sky_mag_arcsec2 = -2.5 * math.log10(sky_brightness)
    # Conversion to nanoLamberts: B_nL = 10^((26.33 - mag) / 2.5)
    sky_nL = 10 ** ((26.33 - sky_mag_arcsec2) / 2.5)

    print("-" * 40)
    print(f"Koefisien Ekstingsi (k)    : {k_zenith:.3f} mag/airmass")
    print(f"Total Ekstingsi @ Bintang  : {ext_corr:.3f} mag")
    print(f"Sky Brightness (mag/arcsec^2): {sky_mag_arcsec2:.2f}")
    print(f"Sky Brightness (nanoLamberts): {sky_nL:.1f} nL")
    print("-" * 40)

    return m_threshold + acuity_corr


if __name__ == "__main__":
    print("Menjalankan Schaefer VISLIMIT High-Precision...")
    
    # 1. Setup Lokasi (Contoh: Semarang)
    lat_semarang = -6.966667
    lon_semarang = 110.416667
    elev_semarang = 10.0
    
    obs = ObservingConditions(
        latitude=lat_semarang,
        longitude=lon_semarang,
        elevation=elev_semarang,
        humidity=70.0,
        temperature=25.0,
        age=23,
        acuity=1.0,
        snellen_factor=6.0
    )

    # 2. Setup Waktu dan Target
    # Tanggal: 16 Januari 2024, 11:30 UT (Sekitar Maghrib di Indonesia)
    yr, mo, dy = 2024, 1, 16 
    ut = 11.5 
    
    # Target (Contoh: Bulan Sabit / Hilal atau Bintang dekat ufuk)
    # Misal star altitude 5 derajat, azimuth 245 derajat
    alt_star = 5.0
    az_star = 245.0

    print(f"\nParameter:")
    print(f"Lokasi: {lat_semarang}, {lon_semarang}")
    print(f"Tanggal: {dy}/{mo}/{yr} jam {ut} UT")
    print(f"Target Alt: {alt_star} deg, Az: {az_star} deg")

    try:
        mag_lim = vislimit_hp(obs, (yr, mo, dy), ut, alt_star, az_star)
        print(f"\nLimiting Magnitude: {mag_lim:.2f}")
    except Exception as e:
        print(f"Error: {e}")