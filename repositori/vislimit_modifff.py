"""
Reconstruction of the VISLIMIT visual limiting magnitude algorithm in Python.

This module provides a programmatic translation of the VISLIMIT calculator
published by Bradley E. Schaefer (Sky & Telescope, May 1998) and later
transcoded into JavaScript by Larry Bogan, with corrections by Victor Reijs.
The implementation here follows the logic of the JavaScript version
available on archaeocosmology.org and replicates its calculations step by
step.  No external astronomy libraries are required; all formulas come
directly from the original open‑source code.

The goal of this module is to expose the same functionality of the web
calculator in a clean Python API.  Two modes are supported:

* **Limiting magnitude mode**: for a given altitude of a target object, the
  algorithm returns the faintest magnitude that can be detected under the
  specified atmospheric conditions.
* **Extinction angle mode**: when a goal magnitude is supplied, the
  algorithm iteratively adjusts the object's topocentric altitude until
  the computed limiting magnitude equals or exceeds the target magnitude.

See the documentation of :func:`vislimit` for details on inputs and
outputs.
"""

from __future__ import annotations

import math
from datetime import date as _date  # for optional date handling in vislimit
from dataclasses import dataclass
from typing import Optional, Tuple


# ---------------------------------------------------------------------------
# Constants for the UBVRI bands
# These arrays store wavelength centres and various constants used in
# extinction and sky‑brightness calculations.  Values are taken directly
# from Schaefer's VISLIMIT program and Bogan's JavaScript version.
WA = (0.365, 0.44, 0.55, 0.70, 0.90)  # microns
MO = (-10.93, -10.45, -11.05, -11.90, -12.70)  # magnitude zero points for Moon
OZ = (0.000, 0.000, 0.031, 0.008, 0.000)  # ozone absorption coefficients
WT = (0.074, 0.045, 0.031, 0.020, 0.015)  # water vapour coefficients
BO = (8.0e-14, 7.0e-14, 1.0e-13, 1.0e-13, 3.0e-13)  # dark sky brightness baseline (ergs/s/cm^2/µm/arcsec^2)
CM = (1.36, 0.91, 0.00, -0.76, -1.17)  # colour correction for Moon brightness
MS = (-25.96, -26.09, -26.74, -27.26, -27.55)  # Sun magnitude in each band


def _deg_to_rad(deg: float) -> float:
    """Convert degrees to radians."""
    return deg * (math.pi / 180.0)


def _rad_to_deg(rad: float) -> float:
    """Convert radians to degrees."""
    return rad * (180.0 / math.pi)


# ---------------------------------------------------------------------------
# DMS<->decimal conversion utilities
#
# Many astronomical inputs (such as latitude and longitude) are often provided
# in degrees, minutes and seconds (DMS).  The VISLIMIT algorithm itself
# operates in decimal degrees internally, but these helper functions make
# it convenient to accept DMS inputs from users.  They mirror the functions
# found in the accompanying Jupyter notebook (repo.ipynb).

def deg_to_des(deg: float, minute: float, second: float) -> float:
    """
    Convert a DMS value to decimal degrees.  The sign of ``deg`` determines
    the sign of the result; ``minute`` and ``second`` should be positive.

    >>> deg_to_des(-8, 40, 0)
    -8.666666666666666

    Args:
        deg: Degrees component (can be negative).
        minute: Minutes component (non‑negative).
        second: Seconds component (non‑negative).

    Returns:
        Decimal degrees with the same sign as ``deg``.
    """
    sign = -1 if deg < 0 else 1
    desimal = abs(deg) + minute / 60.0 + second / 3600.0
    return desimal * sign


def format_dms(desimal: float) -> str:
    """
    Convert decimal degrees to a string in degrees, minutes, and seconds.

    Args:
        desimal: Angle in decimal degrees.

    Returns:
        A string of the form ``"±D° M' S.SS""``.
    """
    sign = -1 if desimal < 0 else 1
    decimal_degrees = abs(desimal)
    degrees = int(decimal_degrees)
    fractional = decimal_degrees - degrees
    total_minutes = fractional * 60.0
    minutes = int(total_minutes)
    seconds = (total_minutes - minutes) * 60.0

    # Apply sign to degrees
    degrees *= sign
    
    return f"{degrees}° {minutes}' {seconds:.2f}\""


def Lintang_tempat(der_lt: float, menit_lt: float, detik_lt: float) -> float:
    """
    Convenience wrapper to convert latitude in DMS to decimal degrees.

    Args:
        der_lt: Degrees of latitude (negative for south).
        menit_lt: Minutes of latitude.
        detik_lt: Seconds of latitude.

    Returns:
        Decimal degrees latitude.
    """
    return deg_to_des(der_lt, menit_lt, detik_lt)


def Bujur_tempat(der_bj: float, menit_bj: float, detik_bj: float) -> float:
    """
    Convenience wrapper to convert longitude in DMS to decimal degrees.

    Args:
        der_bj: Degrees of longitude (negative for west).
        menit_bj: Minutes of longitude.
        detik_bj: Seconds of longitude.

    Returns:
        Decimal degrees longitude.
    """
    return deg_to_des(der_bj, menit_bj, detik_bj)


def _angular_separation(alt1: float, az1: float, alt2: float, az2: float) -> float:
    """
    Compute the angular separation between two points on the celestial sphere
    given their altitudes and azimuths.  All inputs and outputs are in degrees.

    Args:
        alt1: altitude of the first point (degrees)
        az1: azimuth of the first point (degrees, measured eastwards from north)
        alt2: altitude of the second point (degrees)
        az2: azimuth of the second point (degrees, measured eastwards from north)

    Returns:
        Angular separation (degrees) between the two points.
    """
    alt1_rad = _deg_to_rad(alt1)
    alt2_rad = _deg_to_rad(alt2)
    # Note: in the original JavaScript code the azimuth difference is taken
    # directly because azimuths are measured eastwards; we follow the same
    # convention here.
    d_az = _deg_to_rad(az2 - az1)
    # Spherical law of cosines
    cos_sep = (math.sin(alt1_rad) * math.sin(alt2_rad) +
               math.cos(alt1_rad) * math.cos(alt2_rad) * math.cos(d_az))
    # Clamp rounding errors
    cos_sep = max(-1.0, min(1.0, cos_sep))
    return _rad_to_deg(math.acos(cos_sep))


@dataclass
class ObservingConditions:
    """
    Container for the atmospheric and observer parameters used by the
    VISLIMIT algorithm.

    Parameters
    ----------
    humidity : float
        Relative humidity in percent (0–100).
    temperature : float
        Air temperature in degrees Celsius.
    latitude : float
        Observer latitude in decimal degrees (positive north).  Use
        :func:`Lintang_tempat` to convert from DMS if needed.
    longitude : float
        Observer longitude in decimal degrees (positive east).  Although
        longitude is not directly used by the VISLIMIT algorithm, it is
        included here for completeness and future extensions.  Convert
        from DMS with :func:`Bujur_tempat`.
    altitude : float
        Observer elevation above sea level in metres.
    month : int
        Month number (1–12).  Used for seasonal dependence of aerosols.
    year : int
        Year (for dark‑sky seasonal variation).
    snellen_ratio : float, optional
        Ratio of observer visual acuity relative to 20/20.  Defaults to 1.0.
    """

    humidity: float
    temperature: float
    latitude: float
    longitude: float
    altitude: float
    month: int
    year: int
    snellen_ratio: float = 1.0


# -----------------------------
# Core VISLIMIT computation (single altitude)
# -----------------------------

def vislimit(
    conditions: ObservingConditions,
    alt_moon: float,
    az_moon: float,
    alt_sun: float,
    az_sun: float,
    alt_obj: Optional[float] = None,
    goal_magnitude: Optional[float] = None,
    obs_date: Optional["datetime.date"] = None,
    phase_moon_deg: Optional[float] = None,
    step: float = 0.01,
    max_iter: int = 2000,
) -> Tuple[float, float]:
    """
    Compute the limiting magnitude of a star or the extinction angle under
    specified observing conditions.

    This function implements the core of the VISLIMIT algorithm.  In
    addition to the original parameters, it allows the user to specify
    ``obs_date``, a Python :class:`datetime.date` object, from which the
    year, month and approximate lunar phase are derived automatically.
    Alternatively, the lunar phase angle in degrees can be supplied via
    ``phase_moon_deg``.  If both ``obs_date`` and ``phase_moon_deg`` are
    provided, the explicit phase value takes precedence.

    Two modes of operation are supported:

    * **Limiting magnitude mode**: ``goal_magnitude`` is ``None``.  The
      function returns the faintest magnitude visible at a given object
      altitude ``alt_obj`` under the specified conditions.
    * **Extinction angle mode**: ``goal_magnitude`` is a finite value.
      The function iteratively increases the object's altitude until the
      computed limiting magnitude meets or exceeds the goal.  The
      returned altitude indicates when a star of that magnitude becomes
      just visible.

    Args:
        conditions: ``ObservingConditions`` structure with humidity,
            temperature, latitude, longitude, altitude, month, year and
            Snellen ratio.  If ``obs_date`` is given, the date's year
            and month override the ``conditions.month`` and
            ``conditions.year`` values when computing seasonal effects.
        alt_moon: Topocentric altitude of the Moon (degrees).  Positive when
            above the horizon; negative when below.
        az_moon: Topocentric azimuth of the Moon (degrees east of north).
        alt_sun: Topocentric altitude of the Sun (degrees).
        az_sun: Topocentric azimuth of the Sun (degrees east of north).
        alt_obj: Initial altitude of the target object (degrees).  Required
            when computing a limiting magnitude; if omitted in extinction
            angle mode, a starting altitude of 0° is used.
        goal_magnitude: Apparent magnitude of the target object.  If set to
            ``None`` or a very large value (e.g. ``math.inf``), the function
            returns the limiting magnitude for ``alt_obj``.  Otherwise, the
            function returns the altitude at which a star of this
            magnitude is just visible.
        obs_date: Optional :class:`datetime.date` representing the date of
            observation.  If provided, the lunar phase and seasonal
            coefficients are computed from this date.
        phase_moon_deg: Optional explicit phase angle of the Moon in
            degrees.  If supplied, it overrides any phase computed from
            ``obs_date``.
        step: Increment in degrees applied to the object's altitude when
            searching for the extinction angle.  Defaults to 0.01°, matching
            the original JavaScript loop.
        max_iter: Maximum number of iterations when searching for the
            extinction angle.  Prevents infinite loops in pathological cases.

    Returns:
        A tuple ``(limiting_magnitude, altitude)``.  In limiting‑magnitude
        mode, ``altitude`` echoes the input ``alt_obj``.  In
        extinction‑angle mode, ``limiting_magnitude`` equals or exceeds
        ``goal_magnitude``, and ``altitude`` is the computed object altitude.
    """
    # Input sanity checks
    if goal_magnitude is None:
        goal_magnitude = float('inf')

    # If no altitude provided, start from zero altitude (horizon)
    if alt_obj is None:
        alt_obj = 0.0

    # Extract observer conditions
    RH = conditions.humidity
    TE = conditions.temperature
    LA = conditions.latitude
    AL = conditions.altitude
    SN = conditions.snellen_ratio

    # Determine year and month.  If obs_date is provided, use its values;
    # otherwise fall back to the values stored in ObservingConditions.  This
    # allows seasonal and solar‑cycle effects to be computed from a specific
    # date rather than from arbitrary month and year numbers.
    if obs_date is not None:
        Y = obs_date.year
        M = obs_date.month
    else:
        Y = conditions.year
        M = conditions.month

    # Compute the lunar phase angle if not explicitly provided.  When a date
    # is given, use a simple approximation: days since the epoch of a known
    # new moon divided by the synodic period (29.530588853 days) yields the
    # phase fraction.  Multiply by 360° to obtain the phase angle.  This
    # approximation is accurate to within a few degrees for most purposes.
    if phase_moon_deg is None:
        if obs_date is not None:
            # Convert date to Julian Day Number.  We compute at midnight UTC.
            year = obs_date.year
            month = obs_date.month
            day = obs_date.day
            # Julian day calculation based on US Naval Observatory formula
            if month <= 2:
                year -= 1
                month += 12
            A = year // 100
            B = 2 - A + (A // 4)
            jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
            # Epoch of known new moon (2000 Jan 6 18:14 UT) = 2451550.1
            days_since_new = jd - 2451550.1
            phase_fraction = (days_since_new / 29.530588853) % 1.0
            phase_moon_deg = phase_fraction * 360.0
        else:
            # Without a date we cannot compute phase; default to 0° (new moon)
            phase_moon_deg = 0.0

    # RA is a “seasonal angle” used in the aerosol coefficient; the formula
    # uses the month number relative to March and converts to radians.  We
    # convert the difference (month-3)*30° to radians.
    RA = _deg_to_rad((M - 3) * 30.0)
    # SL is the sign of the latitude (±1) to discriminate hemispheres
    SL = 1.0 if LA >= 0 else -1.0
    LT_rad = _deg_to_rad(LA)

    # Limit for the iteration.  We'll loop until we meet the goal magnitude.
    iteration = 0
    limiting_mag = -99.0
    current_alt = alt_obj

    while iteration < max_iter:
        iteration += 1

        # Convert altitudes to zenith distances
        Z = 90.0 - current_alt
        ZM = 90.0 - alt_moon
        ZS = 90.0 - alt_sun

        # Compute angular separations between the object and Moon/Sun
        RM = _angular_separation(current_alt, 0.0, alt_moon, az_moon)  # object azimuth = 0 in this formulation
        RS = _angular_separation(current_alt, 0.0, alt_sun, az_sun)

        # Airmass for each extinction component.  These approximations are
        # empirical fits used in the original program.
        ZZ = _deg_to_rad(Z)
        XG = 1.0 / (math.cos(ZZ) + 0.0286 * math.exp(-10.5 * math.cos(ZZ)))
        XA = 1.0 / (math.cos(ZZ) + 0.0123 * math.exp(-24.5 * math.cos(ZZ)))
        XO_airmass = 1.0 / math.sqrt(1.0 - (math.sin(ZZ) / (1.0 + 20.0 / 6378.0))**2)

        # Containers for extinction and sky brightness per band
        K = [0.0] * 5
        DM = [0.0] * 5
        B = [0.0] * 5

        # Loop over spectral bands UBVRI
        for i in range(5):
            # Extinction coefficients
            KR = 0.1066 * math.exp(-AL / 8200.0) * (WA[i] / 0.55)**(-4.0)
            KA = 0.1 * (WA[i] / 0.55)**(-1.3) * math.exp(-AL / 1500.0)
            # Apply humidity and seasonal corrections.  Victor Reijs
            # introduced a 0.33 factor to the seasonal term; see
            # archaeocosmology.org web page for details【776617556924793†L60-L79】.
            KA *= (1.0 - 0.32 / math.log(RH / 100.0))**1.33 * (1.0 + 0.33 * SL * math.sin(RA))
            KO = OZ[i] * (3.0 + 0.4 * (LT_rad * math.cos(RA) - math.cos(3.0 * LT_rad))) / 3.0
            KW = WT[i] * 0.94 * (RH / 100.0) * math.exp(TE / 15.0) * math.exp(-AL / 8200.0)
            K[i] = KR + KA + KO + KW
            # Total extinction (magnitudes) along the line of sight
            DM[i] = KR * XG + KA * XA + KO * XO_airmass + KW * XG

        # Sky brightness calculations.  Airmasses for star, Moon and Sun
        X = 1.0 / (math.cos(ZZ) + 0.025 * math.exp(-11.0 * math.cos(ZZ)))
        XM = 1.0 / (math.cos(_deg_to_rad(ZM)) + 0.025 * math.exp(-11.0 * math.cos(_deg_to_rad(ZM))))
        if ZM > 90.0:
            XM = 40.0
        XS = 1.0 / (math.cos(_deg_to_rad(ZS)) + 0.025 * math.exp(-11.0 * math.cos(_deg_to_rad(ZS))))
        if ZS > 90.0:
            XS = 40.0

        for i in range(5):
            # Dark night sky baseline with 11‑year solar cycle modulation
            BN = BO[i] * (1.0 + 0.3 * math.cos(6.283 * (Y - 1992) / 11.0))
            BN *= (0.4 + 0.6 / math.sqrt(1.0 - 0.96 * math.sin(ZZ)**2))
            BN *= 10.0 ** (-0.4 * K[i] * X)
            # Moonlight contribution
            MM = -12.73 + 0.026 * abs(phase_moon_deg) + 4e-9 * (phase_moon_deg**4)
            MM += CM[i]
            C3 = 10.0 ** (-0.4 * K[i] * XM)
            # Scattering function FM combines small‑angle (~1/R^2) and wide‑angle terms【776617556924793†L90-L103】
            FM = 6.2e7 / (RM * RM + 1e-10) + 10.0 ** (6.15 - RM / 40.0)
            FM += 10.0 ** 5.36 * (1.06 + math.cos(_deg_to_rad(RM))**2)
            BM = 10.0 ** (-0.4 * (MM - MO[i] + 43.27))
            BM *= (1.0 - 10.0 ** (-0.4 * K[i] * X))
            BM *= (FM * C3 + 440000.0 * (1.0 - C3))
            # Twilight (Sun below horizon) contribution【776617556924793†L104-L106】
            HS = 90.0 - ZS
            BT = 10.0 ** (-0.4 * (MS[i] - MO[i] + 32.5 - HS - (Z / (360.0 * K[i]))))
            BT *= (100.0 / max(RS, 1e-6)) * (1.0 - 10.0 ** (-0.4 * K[i] * X))
            # Daylight contribution【776617556924793†L109-L115】
            C4 = 10.0 ** (-0.4 * K[i] * XS)
            FS = 6.2e7 / (RS * RS + 1e-10) + 10.0 ** (6.15 - RS / 40.0)
            FS += 10.0 ** 5.36 * (1.06 + math.cos(_deg_to_rad(RS))**2)
            BD = 10.0 ** (-0.4 * (MS[i] - MO[i] + 43.27))
            BD *= (1.0 - 10.0 ** (-0.4 * K[i] * X))
            BD *= (FS * C4 + 440000.0 * (1.0 - C4))
            # Choose the smaller of twilight and daylight when both Sun and Moon below horizon
            sky_component = BD if BD < BT else BT
            total = BN + sky_component
            # Add moonlight if Moon is above horizon
            if ZM < 90.0:
                total += BM
            # Convert to picoergs
            B[i] = total * 1e12

        # Convert band V brightness to nanolamberts; 1 nL = 1.02e-3 picoergs/s/cm^2/µm/arcsec^2
        BL = B[2] / 1.02e-3
        # Select constants based on brightness threshold【776617556924793†L124-L134】
        if BL < 1500.0:
            C1 = 10.0 ** (-9.8)
            C2 = 10.0 ** (-1.9)
        else:
            C1 = 10.0 ** (-8.350001)
            C2 = 10.0 ** (-5.9)
        TH = C1 * (1.0 + math.sqrt(C2 * BL))**2
        # Compute limiting magnitude with Snellen ratio correction【776617556924793†L134-L136】
        limiting_mag = -16.57 - 2.5 * math.log10(TH) - DM[2] + 5.0 * math.log10(SN)

        # Check goal/iteration conditions
        if limiting_mag >= goal_magnitude - 1e-6:
            # Either returning after computing magnitude (goal is infinite) or reaching target magnitude
            return limiting_mag, current_alt

        # Not yet at goal magnitude: increment altitude and repeat
        current_alt += step

    # If we exit the loop without meeting the goal, return the last computed values
    return limiting_mag, current_alt


def example_usage():
    """
    Example invocation of the vislimit algorithm.  This demonstrates how to
    compute the limiting magnitude for a star 10° above the horizon given
    certain observing conditions, or determine the altitude at which a star
    of magnitude 4.0 becomes just visible.

    The numbers in this function are arbitrary; users should replace them
    with real observational parameters.
    """
    # Define observing conditions
    cond = ObservingConditions(
        humidity=92,
        temperature=24,
        latitude=-10,
        longitude=0.0,
        altitude=3,
        month=1,
        year=2024,
        snellen_ratio=1.0,
    )

    # Compute limiting magnitude for a star at 19° altitude
    from datetime import date

    lm, alt = vislimit(
        conditions=cond,
        alt_moon=-90.0,
        az_moon=180,
        alt_sun=-0.3,
        az_sun=0.142,
        alt_obj=11.83,
        goal_magnitude=99,
        obs_date=date(2024, 1, 11),
    )
    print(f"Limiting magnitude at {format_dms(19.05)}: {lm:.2f}")

    # Determine altitude at which a magnitude‑4 star becomes visible
    lm_goal, alt_goal = vislimit(
        conditions=cond,
        alt_moon=-90.0,
        az_moon=0.0,
        alt_sun=0.0,
        az_sun=0.0,
        alt_obj=0.0,
        goal_magnitude=4.0,
        obs_date=date(1972, 3, 15),
    )
    print(f"A magnitude‑4 star becomes visible at {format_dms(alt_goal)} altitude; limiting mag = {lm_goal:.2f}")


if __name__ == "__main__":
    example_usage()