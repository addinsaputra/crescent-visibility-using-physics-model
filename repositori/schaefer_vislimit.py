"""
Implementation of Bradley E. Schaefer’s visual limiting magnitude algorithm in
Python.

This module provides a function, ``vislimit``, which computes the limiting
visual magnitude for a star given the observing conditions.  The algorithm is
adapted from the BASIC program VISLIMIT.BAS published in _Sky & Telescope_,
May 1998, and from the JavaScript translation by Larry Bogan.  The core
calculations follow the formulas documented in the original code.  For
references to the specific steps see the cited lines in Bogan’s JavaScript
version【665925694870024†L309-L315】【665925694870024†L339-L356】.

Key features:

* Conversion between calendar date and Julian Date.
* Local Sidereal Time computation following the Astronomical Almanac
  convention【665925694870024†L267-L299】.
* Low‑precision Sun and Moon positions with topocentric correction
  taken from the Astronomical Almanac【665925694870024†L309-L315】【665925694870024†L339-L356】.
* Airmass, extinction, and sky brightness calculations for the V band
  using the formulas and constants from Schaefer’s program【665925694870024†L455-L486】.
* Phase function and glare corrections for the Moon and Sun【665925694870024†L507-L515】.
* Twilight, night and daylight sky brightness contributions【665925694870024†L561-L566】【665925694870024†L571-L575】【665925694870024†L579-L585】.
* Pupil size, critical visual angle and threshold of visibility adjustments
  based on observer age and atmospheric conditions【665925694870024†L640-L667】.

The function returns the limiting magnitude in the V band.  Users can
instantiate their own observing conditions and call ``vislimit`` repeatedly
to build tables of visibility for different dates and times.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def julian_date(year: int, month: int, day: int, hour: float) -> float:
    """
    Compute the Julian Date for a calendar date and fractional hour.

    The algorithm below follows the standard astronomical formula and does not
    apply any delta‑T corrections.  The input ``hour`` should be the UTC time
    expressed in decimal hours.  For example, 18:30 UTC would be ``18.5``.

    Args:
        year: Gregorian year (positive integer).
        month: 1–12.
        day: day of the month.
        hour: fractional hour (UTC).

    Returns:
        The Julian Date (days since 1 Jan 4713 BC noon UTC).
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
    """
    Compute the local mean sidereal time (in hours) at the given Julian Date
    and east longitude.【665925694870024†L267-L299】

    Args:
        jd: Julian Date.
        longitude_deg: East longitude in degrees (positive eastwards, negative west).

    Returns:
        Sidereal time in decimal hours between 0 and 24.
    """
    # Extract integer and fractional parts of JD to get UT
    jd_int = math.floor(jd)
    jd_frac = jd - jd_int
    if jd_frac < 0.5:
        jd_mid = jd_int - 0.5
        ut = jd_frac + 0.5
    else:
        jd_mid = jd_int + 0.5
        ut = jd_frac - 0.5
    t = (jd_mid - 2451545.0) / 36525.0  # J2000 = 2451545.0
    # Sidereal time at Greenwich (fraction of a day)
    sid_g = (24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 6.2e-6 * t * t * t) / 86400.0
    sid_g -= math.floor(sid_g)  # fractional part
    sid_g += 1.0027379093 * ut - longitude_deg / 360.0
    sid_g = (sid_g - math.floor(sid_g)) * 24.0
    if sid_g < 0:
        sid_g += 24.0
    return sid_g


def atan_circ(x: float, y: float) -> float:
    """
    Compute the full‑circle arctangent of y/x returning an angle in the range
    [0, 2π).【665925694870024†L69-L82】

    Args:
        x: x‑coordinate.
        y: y‑coordinate.

    Returns:
        Angle in radians measured from the positive x‑axis.
    """
    if x == 0.0:
        if y > 0.0:
            return math.pi / 2
        elif y < 0.0:
            return 3 * math.pi / 2
        else:
            return 0.0
    angle = math.atan(y / x)
    if x < 0.0:
        angle += math.pi
    if angle < 0.0:
        angle += 2 * math.pi
    return angle


def low_precision_sun(jd: float) -> Tuple[float, float]:
    """
    Compute the Sun’s right ascension (hours) and declination (degrees) using
    low‑precision formulas.【665925694870024†L309-L322】

    Args:
        jd: Julian Date.

    Returns:
        (ra_sun_hours, dec_sun_degrees)
    """
    n = jd - 2451545.0  # days since J2000
    L = 280.460 + 0.9856474 * n
    g = (357.528 + 0.9856003 * n) * (math.pi / 180.0)
    lambda_sun = (L + 1.915 * math.sin(g) + 0.020 * math.sin(2 * g)) * (math.pi / 180.0)
    epsilon = (23.439 - 0.0000004 * n) * (math.pi / 180.0)
    x = math.cos(lambda_sun)
    y = math.cos(epsilon) * math.sin(lambda_sun)
    z = math.sin(epsilon) * math.sin(lambda_sun)
    ra = atan_circ(x, y) * (12.0 / math.pi)  # radians to hours (multiply by 12/π)
    dec = math.asin(z) * (180.0 / math.pi)
    return ra, dec


def low_precision_moon(jd: float, lat_deg: float, lst_hours: float) -> Tuple[float, float, float]:
    """
    Compute the Moon’s topocentric right ascension (hours), declination (degrees)
    and distance (Earth radii) using low‑precision series.【665925694870024†L339-L356】

    Args:
        jd: Julian Date.
        lat_deg: Observer latitude in degrees.
        lst_hours: Local sidereal time in hours.

    Returns:
        (ra_moon_hours, dec_moon_degrees, dist_moon_earth_radii)
    """
    T = (jd - 2451545.0) / 36525.0
    # Ecliptic longitude and latitude (degrees)
    lambda_moon = (
        218.32
        + 481267.883 * T
        + 6.29 * math.sin((134.9 + 477198.85 * T) * math.pi / 180.0)
        - 1.27 * math.sin((259.2 - 413335.38 * T) * math.pi / 180.0)
        + 0.66 * math.sin((235.7 + 890534.23 * T) * math.pi / 180.0)
        + 0.21 * math.sin((269.9 + 954397.70 * T) * math.pi / 180.0)
        - 0.19 * math.sin((357.5 + 35999.05 * T) * math.pi / 180.0)
        - 0.11 * math.sin((186.6 + 966404.05 * T) * math.pi / 180.0)
    )
    beta = (
        5.13 * math.sin((93.3 + 483202.03 * T) * math.pi / 180.0)
        + 0.28 * math.sin((228.2 + 960400.87 * T) * math.pi / 180.0)
        - 0.28 * math.sin((318.3 + 6003.18 * T) * math.pi / 180.0)
        - 0.17 * math.sin((217.6 - 407332.20 * T) * math.pi / 180.0)
    )
    # Convert to radians
    lambda_rad = lambda_moon * math.pi / 180.0
    beta_rad = beta * math.pi / 180.0
    # Parallax (pie) and distance
    pie = (
        0.9508
        + 0.0518 * math.cos((134.9 + 477198.85 * T) * math.pi / 180.0)
        + 0.0095 * math.cos((259.2 - 413335.38 * T) * math.pi / 180.0)
        + 0.0078 * math.cos((235.7 + 890534.23 * T) * math.pi / 180.0)
        + 0.0028 * math.cos((269.9 + 954397.70 * T) * math.pi / 180.0)
    )
    pie_rad = pie * math.pi / 180.0
    dist = 1.0 / math.sin(pie_rad)
    # Geocentric rectangular coordinates
    l = math.cos(beta_rad) * math.cos(lambda_rad)
    m = 0.9175 * math.cos(beta_rad) * math.sin(lambda_rad) - 0.3978 * math.sin(beta_rad)
    n = 0.3978 * math.cos(beta_rad) * math.sin(lambda_rad) + 0.9175 * math.sin(beta_rad)
    x = l * dist
    y = m * dist
    z = n * dist
    # Topocentric correction
    lat_rad = lat_deg * math.pi / 180.0
    lst_rad = lst_hours * (math.pi / 12.0)  # hours to radians (24h->2π)
    x -= math.cos(lat_rad) * math.cos(lst_rad)
    y -= math.cos(lat_rad) * math.sin(lst_rad)
    z -= math.sin(lat_rad)
    topo_dist = math.sqrt(x * x + y * y + z * z)
    l = x / topo_dist
    m = y / topo_dist
    n = z / topo_dist
    alpha = atan_circ(l, m)
    delta = math.asin(n)
    ra = alpha * (12.0 / math.pi)  # convert to hours
    dec = delta * (180.0 / math.pi)
    return ra, dec, topo_dist


def altaz_to_radec(
    alt_deg: float, az_deg: float, lat_deg: float, lst_hours: float
) -> Tuple[float, float]:
    """
    Convert an altitude/azimuth pair to equatorial coordinates (RA in hours,
    Dec in degrees) for a given latitude and local sidereal time.  The
    azimuth is measured eastwards from north.

    Args:
        alt_deg: Altitude in degrees above the horizon.
        az_deg: Azimuth in degrees east of north.
        lat_deg: Observer latitude in degrees.
        lst_hours: Local sidereal time in hours.

    Returns:
        (ra_hours, dec_degrees).
    """
    alt_rad = math.radians(alt_deg)
    az_rad = math.radians(az_deg)
    lat_rad = math.radians(lat_deg)
    # Compute declination
    sin_dec = math.sin(alt_rad) * math.sin(lat_rad) + math.cos(alt_rad) * math.cos(lat_rad) * math.cos(az_rad)
    dec_rad = math.asin(sin_dec)
    # Hour angle via azimuth formula
    # Using formula: tan(H) = sin(A) / (cos(A)*sin(lat) + tan(alt)*cos(lat))
    ha_rad = math.atan2(
        -math.sin(az_rad) * math.cos(alt_rad),
        math.cos(az_rad) * math.cos(alt_rad) * math.sin(lat_rad) - math.tan(alt_rad) * math.cos(lat_rad),
    )
    # Convert to RA
    ha_hours = ha_rad * (12.0 / math.pi)
    ra_hours = lst_hours - ha_hours
    # Normalize RA to [0,24)
    ra_hours = ra_hours % 24.0
    dec_deg = math.degrees(dec_rad)
    return ra_hours, dec_deg


@dataclass
class ObservingConditions:
    """Container for observer and environmental parameters."""

    latitude: float  # degrees
    longitude: float  # degrees (east positive)
    elevation: float  # metres above sea level
    humidity: float  # relative humidity in percent (0–100)
    temperature: float  # degrees Celsius
    age: float  # observer age in years
    acuity: float  # visual acuity (Snellen fraction; 6 means normal)
    snellen_factor: float  # Snellen n/r (value used as qsneln)


def vislimit(
    obs: ObservingConditions,
    date: Tuple[int, int, int],
    ut_hour: float,
    alt_star: float,
    az_star: float,
) -> float:
    """
    Compute the visual limiting magnitude for a star at the given altitude and
    azimuth under the specified observing conditions.

    The calculation follows Schaefer’s VISLIMIT algorithm.  Only the V band
    (wavelength 0.55 µm) is considered, because this band best matches the
    human eye’s photopic response【665925694870024†L455-L486】.

    Args:
        obs: ObservingConditions containing site and observer parameters.
        date: Tuple of (year, month, day) in the Gregorian calendar.
        ut_hour: Universal Time as a decimal hour (e.g. 13.5 for 13:30 UT).
        alt_star: Altitude of the target star in degrees (0–90).  Values outside
            this range will be clipped.
        az_star: Azimuth of the star in degrees east of north (0–360).

    Returns:
        The limiting visual magnitude in the V band.
    """
    year, month, day = date
    # Clip star altitude to [0, 90]
    alt_star = max(0.0, min(90.0, alt_star))
    # Convert to Julian Date
    jd = julian_date(year, month, day, ut_hour)
    # Compute local sidereal time
    lst = local_sidereal_time(jd, obs.longitude)
    # Solar position
    ra_sun, dec_sun = low_precision_sun(jd)
    # Lunar position
    ra_moon, dec_moon, dist_moon = low_precision_moon(jd, obs.latitude, lst)
    # Convert star alt/az to RA/Dec
    ra_star, dec_star = altaz_to_radec(alt_star, az_star, obs.latitude, lst)
    # Convert altitudes for Sun and Moon
    # alt = arcsin(sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(HA))
    def equ_to_alt(ra_h: float, dec_d: float) -> float:
        ha_h = (lst - ra_h + 24.0) % 24.0
        ha_rad = ha_h * (math.pi / 12.0)
        dec_rad = math.radians(dec_d)
        lat_rad = math.radians(obs.latitude)
        alt_rad = math.asin(
            math.sin(dec_rad) * math.sin(lat_rad)
            + math.cos(dec_rad) * math.cos(lat_rad) * math.cos(ha_rad)
        )
        return math.degrees(alt_rad)
    alt_sun = equ_to_alt(ra_sun, dec_sun)
    alt_moon = equ_to_alt(ra_moon, dec_moon)
    # Zenith distances (radians)
    zz_star = math.radians(90.0 - alt_star)
    zz_moon = math.radians(90.0 - alt_moon)
    zz_sun = math.radians(90.0 - alt_sun)
    # Airmass for extinction (star)【665925694870024†L455-L486】
    xg = 1.0 / (math.cos(zz_star) + 0.0286 * math.exp(-10.5 * math.cos(zz_star)))
    xa = 1.0 / (math.cos(zz_star) + 0.0123 * math.exp(-24.5 * math.cos(zz_star)))
    xo = 1.0 / math.sqrt(1.0 - (math.sin(zz_star) / (1.0 + 20.0 / 6378.0)) ** 2)
    # Wavelength index (3rd element, index 2) for V band
    wa = 0.55  # microns
    oz = 0.031  # ozone coefficient for V band【665925694870024†L397-L398】
    wt = 0.031  # water vapour coefficient for V band【665925694870024†L398-L399】
    # Extinction coefficients【665925694870024†L475-L483】
    kr = 0.1066 * math.exp(-obs.elevation / 8200.0) * (wa / 0.55) ** (-4.0)
    # Avoid log of zero; humidity as fraction
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
    kw = (
        wt
        * 0.94
        * humid_frac
        * math.exp(obs.temperature / 15.0)
        * math.exp(-obs.elevation / 8200.0)
    )
    k_tot = kr + ka + ko + kw
    # Differential magnitude due to extinction【665925694870024†L475-L486】
    dm_tot = kr * xg + ka * xa + ko * xo + kw * xg
    # Airmass for Sun and Moon【665925694870024†L494-L502】
    def airmass(zr: float) -> float:
        am = 1.0 / (math.cos(zr) + 0.025 * math.exp(-11.0 * math.cos(zr)))
        # If object is below horizon (z > 90°), set airmass to large value
        if zr > math.pi / 2:
            return 40.0
        return am
    am_star = airmass(zz_star)
    am_moon = airmass(zz_moon)
    am_sun = airmass(zz_sun)
    # Angular separations (radians)【665925694870024†L431-L435】
    def subtend(ra1_h: float, dec1_d: float, ra2_h: float, dec2_d: float) -> float:
        ra1 = ra1_h * (math.pi / 12.0)
        dec1 = math.radians(dec1_d)
        ra2 = ra2_h * (math.pi / 12.0)
        dec2 = math.radians(dec2_d)
        # Use spherical law of cosines
        cos_theta = (
            math.sin(dec1) * math.sin(dec2)
            + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
        )
        # Clamp due to numerical errors
        cos_theta = max(-1.0, min(1.0, cos_theta))
        theta = math.acos(cos_theta)
        return theta
    # Angular distances
    rho_moon = subtend(ra_star, dec_star, ra_moon, dec_moon)
    rho_sun = subtend(ra_star, dec_star, ra_sun, dec_sun)
    rho_sm = subtend(ra_moon, dec_moon, ra_sun, dec_sun)
    # Illumination fraction of the Moon (percentage).  The phase angle is the
    # Sun–Moon separation.  Illumination fraction f = (1 + cos(phase))/2.
    f_illum = 50.0 * (1.0 + math.cos(rho_sm))
    # Phase angle qphase (in radians)【665925694870024†L429-L434】
    qphase = math.acos(2.0 * f_illum / 100.0 - 1.0)
    # Phase functions for Moon and Sun【665925694870024†L507-L515】
    def phase_function(rho: float) -> float:
        return (
            10.0 ** 5.36 * (1.06 + math.cos(rho) ** 2)
            + 10.0 ** (6.15 - rho * (180.0 / math.pi) / 40.0)
            + 6.2e7 * (rho * (180.0 / math.pi)) ** (-2.0)
        )
    ff_moon = phase_function(rho_moon)
    ff_sun = phase_function(rho_sun)
    ff_sunm = phase_function(rho_sm)
    # Moon brightness【665925694870024†L519-L531】
    zmag_mn = -12.73 + 0.026 * abs(qphase) + 4e-9 * (abs(qphase) ** 4)
    is_moon = 10.0 ** (-0.4 * (zmag_mn + 16.57))
    # Distance scaling (mean distance 60.27 Earth radii)
    is_moon /= (dist_moon / 60.27) ** 2
    # Opposition effect
    is_moon *= max(1.0, (1.35 - 0.05 * abs(qphase)))
    # Apply extinction to moon brightness at lunar airmass
    i_moon = is_moon * 10.0 ** (-0.4 * k_tot * am_moon)
    # Surface brightness and glare【665925694870024†L542-L556】
    if rho_moon * (180.0 / math.pi) <= 5.0 and qphase > 1.0:
        # Moon glare for small angular separations
        b_atm = 6.25e7 * is_moon * (rho_moon) ** (-2.0) * (
            10.0 ** (-0.4 * k_tot * am_moon)
            - 10.0 ** (-0.8 * k_tot * am_moon)
        )
        b_eye = 4.63e7 * i_moon * (rho_moon) ** (-2.0)
        b_glare = b_atm + b_eye
        b_moon = 5.67e10 * i_moon / qphase
    else:
        b_glare = 0.0
        b_moon = ff_moon * i_moon * (1.0 - 10.0 ** (-0.4 * k_tot * am_star))
    # Night sky brightness【665925694870024†L561-L566】
    # bo is 200 for V band (index 3)【665925694870024†L399-L401】
    bo_v = 200.0
    year = date[0]
    b_night = bo_v * (0.4 + 0.6 * am_star) * 10.0 ** (-0.4 * k_tot * am_star) * (
        1.0 + 0.3 * math.cos(2.0 * math.pi * (year - 1992.0) / 11.0)
    )
    b_night_m = bo_v * (0.4 + 0.6 * am_moon) * 10.0 ** (-0.4 * k_tot * am_moon) * (
        1.0 + 0.3 * math.cos(2.0 * math.pi * (year - 1992.0) / 11.0)
    )
    # Twilight brightness【665925694870024†L571-L575】
    h_sun = alt_sun
    # Ensure the exponent does not overflow
    # Twilight term: max(1, 10^(rho_sun*deg/90 - 1.1)) * 10^(8.45 + 0.4*h_sun)
    btwilt = max(1.0, 10.0 ** (rho_sun * (180.0 / math.pi) / 90.0 - 1.1)) * 10.0 ** (
        8.45 + 0.4 * h_sun
    ) * (1.0 - 10.0 ** (-0.4 * k_tot * am_star))
    btwilt_m = max(1.0, 10.0 ** (rho_sm * (180.0 / math.pi) / 90.0 - 1.1)) * 10.0 ** (
        8.45 + 0.4 * h_sun
    ) * (1.0 - 10.0 ** (-0.4 * k_tot * am_moon))
    # Daylight brightness【665925694870024†L579-L585】
    b_day = (
        11700.0
        * ff_sun
        * 10.0 ** (-0.4 * k_tot * am_sun)
        * (1.0 - 10.0 ** (-0.4 * k_tot * am_star))
    )
    b_day_m = (
        11700.0
        * ff_sunm
        * 10.0 ** (-0.4 * k_tot * am_sun)
        * (1.0 - 10.0 ** (-0.4 * k_tot * am_moon))
    )
    # Effective background brightness【665925694870024†L598-L604】
    beff = b_moon + b_glare + b_night + min(btwilt, b_day)
    beff_m = b_night_m + min(btwilt_m, b_day_m)
    # Day/night classification (above or below threshold 1479 nL)【665925694870024†L610-L617】
    day = 1 if beff >= 1479.0 else 0
    day_m = 1 if beff_m >= 1479.0 else 0
    # Correction factors【665925694870024†L620-L633】
    if day == 1:
        Fe = 10.0 ** (0.4 * k_tot * am_star)
        Fci = 1.0
        Fcb = 1.0
    else:
        Fe = 10.0 ** (0.48 * k_tot * am_star)
        Fci = 10.0 ** (-0.4 * (1.0 - 0.5 / 2.0))
        Fcb = 10.0 ** (-0.4 * (1.0 - 0.7 / 2.0))
    # Pupil size factor【665925694870024†L640-L645】
    De = 7.0 * math.exp(-0.5 * (obs.age / 100.0) ** 2)
    Ds = 7.0 * math.exp(-0.5 * (25.0 / 100.0) ** 2)
    Fp = (Ds / De) ** 2
    # Critical visual angle【665925694870024†L648-L667】
    qsneln = max(1e-6, obs.snellen_factor)
    if day == 1:
        cva = 42.0 * 10.0 ** (8.28 * (beff ** -0.29)) / qsneln
    else:
        cva = min(900.0, 380.0 * 10.0 ** (0.3 * (beff ** -0.29))) / qsneln
    zeta = 1800.0 * math.sqrt(f_illum)  # equivalent angular diameter of moon
    if cva > zeta:
        cth = max(2.4 * (beff ** -0.1), 20.0 * (beff ** -0.4)) * ((cva / zeta) ** 2)
    elif beff > 1.0e6:
        cth = 0.0028 + 2.4 * (beff ** -0.1) * ((cva / zeta) ** 2)
    else:
        cth = 10.0 ** (
            -1.0 * (0.12 * 0.40 * math.log(beff) / math.log(10.0))
            + (0.90 - 0.15 * math.log(beff) / math.log(10.0)) * math.log(100.0 * 60.0 / zeta) / math.log(10.0)
        )
    # Seeing correction【665925694870024†L676-L681】
    xi = math.sqrt(8.0 * 0.361 * (1.5 * 1.5 * am_star))
    Fr = (1.0 + 0.03 * (xi / cva) ** 2) / (qsneln ** 2)
    # Effective brightness perceived by the eye【665925694870024†L735-L748】
    b_effective = beff / (Fp * Fcb)
    if day == 1:
        c1 = 10.0 ** -8.35
        c2 = 10.0 ** -5.9
    else:
        c1 = 10.0 ** -9.8
        c2 = 10.0 ** -1.9
    zlim_fc = c1 * (1.0 + math.sqrt(c2 * b_effective)) ** 2
    # Illuminance outside the atmosphere
    zlim_fcs = zlim_fc * (Fp * Fr * Fci * Fe)
    # Convert to magnitudes【665925694870024†L746-L759】
    zlim_mag = -16.57 - 2.5 * math.log10(zlim_fcs)
    # Apply extinction and acuity corrections【665925694870024†L759-L763】
    zlim = zlim_mag + 0.16 * (obs.acuity - 6.0)
    return zlim - dm_tot