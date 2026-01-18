"""VISLIMIT (Schaefer/Bogan/Reijs) reconstruction in Python.

Goal
----
Match the *outputs* and *workflow* of the web calculator:
https://www.archaeocosmology.org/eng/vislimit-Bogan-vr.htm

This module reconstructs the Schaefer (1998) visual limiting magnitude
algorithm as implemented by Larry Bogan (JavaScript) and corrected/extended
by Victor Reijs.

Important conventions (same as the Reijs web page)
--------------------------------------------------
- Inputs for Moon/Sun azimuth are **"Azimuth compared to Object"** (degrees).
  That means the object is treated as azimuth=0 and the Moon/Sun are at the
  given relative azimuth.

- "Goal Magnitude" controls the mode:
  * Goal Magnitude = 99 (or None in this Python API) => compute limiting
    magnitude at the given object altitude.
  * Goal Magnitude != 99 => compute the object altitude where the object of
    that magnitude becomes just visible (altitude search).

Units
-----
- Angles are in **degrees**.
- Altitudes are topocentric altitudes.
- Sky brightness is reported as:
  * V-band nanoLamberts (nL) in the RESULTS section.
  * UBVRI picoergs/(cm^2*micron*arcsec^2) in Supplemental Results.

Modern date/phase option
------------------------
The original web page asks for Month/Year + a coarse "Phase of the Moon"
(New/Crescent/Quarter/Gibbous/Full). For scientific use you usually want a
continuous lunar phase angle.

This module accepts either:
- phase_moon_deg: explicit phase angle in degrees, or
- obs_datetime_utc: a datetime (UTC) from which we estimate lunar phase angle.

If Skyfield is installed, the phase angle is computed from ephemeris
(sun-moon elongation as seen from Earth). If not, we fall back to a standard
synodic-period approximation.

References
----------
- Bradley E. Schaefer, "To the Visual Limits", Sky & Telescope, May 1998.
- Reijs web page and log of JS typos:
  https://www.archaeocosmology.org/eng/vislimit-Bogan-vr.htm
  https://www.archaeocosmology.org/eng/vislimitLogan.htm
"""

from __future__ import annotations #For forward compatibility

import math #Math functions
from dataclasses import dataclass #Data classes
from datetime import date, datetime, timezone #Date and time
from typing import Any, Dict, Optional, Tuple #Type hints

# -----------------------------
# Constants used in Schaefer/Bogan/Reijs
# -----------------------------

# Effective wavelengths (microns)
WAC = (0.36, 0.44, 0.55, 0.64, 0.79)

# Ozone absorption coefficients
OZ = (0.0, 0.0, 0.031, 0.008, 0.0)

# Water vapor absorption coefficients
WT = (0.074, 0.045, 0.031, 0.020, 0.015)

# Baseline dark sky brightness (erg/s/cm^2/µm/arcsec^2)
BO = (6.05e-14, 2.56e-13, 1.02e-13, 8.10e-14, 1.23e-13)

# Band offsets for Moon magnitude
CM = (0.0, 0.0, 0.0, -0.10, -0.30)

# Magnitude zero points for Moon
MO = (-10.93, -10.45, -11.05, -11.90, -12.70)

# Sun magnitudes per band
MS = (-25.96, -26.09, -26.74, -27.26, -27.55)

_BANDS = ("U", "B", "V", "R", "I")


@dataclass(frozen=True) #Immutable data class
class ObservingConditions:
    humidity_percent: float #Relative humidity
    temperature_c: float #Temperature in Celsius
    latitude_deg: float #Latitude in degrees
    height_m: float #Altitude in meters
    snellen_ratio: float = 1.0 #Snellen ratio (20/20 => 1; 20/10 => 2)


# -----------------------------
# Helpers
# -----------------------------

def _deg_to_rad(x_deg: float) -> float:
    return x_deg * math.pi / 180.0


def _angular_separation_deg(alt1_deg: float, az1_deg: float, alt2_deg: float, az2_deg: float) -> float:
    """Great-circle separation on the unit sphere from alt/az (degrees)."""
    alt1 = _deg_to_rad(alt1_deg)
    alt2 = _deg_to_rad(alt2_deg)
    daz = _deg_to_rad(az1_deg - az2_deg)
    cos_sep = math.sin(alt1) * math.sin(alt2) + math.cos(alt1) * math.cos(alt2) * math.cos(daz)
    cos_sep = max(-1.0, min(1.0, cos_sep))
    return math.degrees(math.acos(cos_sep))


def _julian_day_utc0(d: date) -> float:
    """Julian Day at 0h UTC for a Gregorian date."""
    y, m, dd = d.year, d.month, d.day
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + (A // 4)
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + dd + B - 1524.5
    return float(jd)


def _moon_phase_angle_deg_from_datetime(obs_datetime_utc: datetime) -> float:
    """Compute a *continuous* lunar phase angle in degrees.

    Uses Skyfield if available, otherwise falls back to a synodic approximation.

    Returned angle is in [0, 360): 0=new, 180=full.
    """
    if obs_datetime_utc.tzinfo is None:
        obs_datetime_utc = obs_datetime_utc.replace(tzinfo=timezone.utc)
    else:
        obs_datetime_utc = obs_datetime_utc.astimezone(timezone.utc)

    # Try Skyfield (preferred)
    try:
        from skyfield.api import load

        ts = load.timescale()
        t = ts.from_datetime(obs_datetime_utc)
        eph = load("de421.bsp")
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]
        e = earth.at(t)
        # Elongation Sun-Moon as seen from Earth in degrees
        elong = e.observe(sun).separation_from(e.observe(moon)).degrees
        # Convert elongation to phase angle in [0,360)
        # elongation is in [0,180]. We need waxing/waning; approximate using ecliptic longitudes.
        # Compute ecliptic longitudes and take signed difference.
        sun_lon, _, _ = e.observe(sun).ecliptic_latlon()
        moon_lon, _, _ = e.observe(moon).ecliptic_latlon()
        dlon = (moon_lon.degrees - sun_lon.degrees) % 360.0
        return float(dlon)
    except Exception:
        # Fallback: synodic approximation using a standard epoch new moon
        # Epoch: 2000-01-06 18:14 UTC ~ JD 2451550.1 (commonly used in simple phase formulas)
        jd = _julian_day_utc0(obs_datetime_utc.date())
        days_since = jd - 2451550.1
        phase_fraction = (days_since / 29.530588853) % 1.0
        return float(phase_fraction * 360.0)


# -----------------------------
# Core VISLIMIT computation (single altitude)
# -----------------------------

def _compute_vislimit_once(
    *,
    cond: ObservingConditions, #Observing conditions
    month: int, #Month
    year: int, #Year
    phase_moon_deg: float, #Moon phase angle in degrees
    alt_obj_deg: float, #Object altitude in degrees
    alt_moon_deg: float,
    az_moon_rel_deg: float,
    alt_sun_deg: float,
    az_sun_rel_deg: float,
) -> Dict[str, Any]:
    """Compute all intermediate quantities + limiting magnitude at one altitude."""

    RH = cond.humidity_percent #relative humidity
    TE = cond.temperature_c #temperature in Celsius
    LA = cond.latitude_deg #latitude in degrees
    AL = cond.height_m #altitude in meters
    SN = cond.snellen_ratio #Snellen ratio

    # Seasonal angle (month relative to March)
    RA = _deg_to_rad((month - 3) * 30.0) #seasonal angle in radians
    SL = 1.0 if LA >= 0 else -1.0 #sign of latitude
    LT = _deg_to_rad(LA) #latitude in radians

    # Geometry
    Z = 90.0 - alt_obj_deg
    ZM = 90.0 - alt_moon_deg
    ZS = 90.0 - alt_sun_deg

    # Object is at az=0; Moon/Sun are given as "azimuth compared to object"
    RM = _angular_separation_deg(alt_obj_deg, 0.0, alt_moon_deg, az_moon_rel_deg)
    RS = _angular_separation_deg(alt_obj_deg, 0.0, alt_sun_deg, az_sun_rel_deg)

    ZZ = _deg_to_rad(Z)

    # Airmass components (as per Bogan/Reijs)
    XG = 1.0 / (math.cos(ZZ) + 0.0286 * math.exp(-10.5 * math.cos(ZZ)))
    XA = 1.0 / (math.cos(ZZ) + 0.0123 * math.exp(-24.5 * math.cos(ZZ)))
    XO = 1.0 / math.sqrt(1.0 - (math.sin(ZZ) / (1.0 + 20.0 / 6378.0)) ** 2)

    # Extinction coefficients and extinction magnitudes
    K: Dict[str, float] = {}
    DM: Dict[str, float] = {}

    for v, band in enumerate(_BANDS):
        # Rayleigh
        KR = 0.1066 * math.exp(-AL / 8200.0) * (WAC[v] / 0.55) ** (-4.0)
        # Aerosol
        KA = 0.1 * (WAC[v] / 0.55) ** (-1.3) * math.exp(-AL / 1500.0)
        KA *= (1.0 - 0.32 / math.log(RH / 100.0)) ** 1.33 * (1.0 + 0.33 * SL * math.sin(RA))
        # Ozone
        KO = OZ[v] * (3.0 + 0.4 * (LT * math.cos(RA) - math.cos(3.0 * LT))) / 3.0
        # Water vapor
        KW = WT[v] * 0.94 * (RH / 100.0) * math.exp(TE / 15.0) * math.exp(-AL / 8200.0)

        k = KR + KA + KO + KW
        K[band] = k
        DM[band] = KR * XG + KA * XA + KO * XO + KW * XG

    # Airmass for sky-brightness scattering model
    X = 1.0 / (math.cos(ZZ) + 0.025 * math.exp(-11.0 * math.cos(ZZ)))

    def _airmass_for_zenith_distance(z_deg: float) -> float:
        z = _deg_to_rad(z_deg)
        xx = 1.0 / (math.cos(z) + 0.025 * math.exp(-11.0 * math.cos(z)))
        return 40.0 if z_deg > 90.0 else xx

    XM = _airmass_for_zenith_distance(ZM)
    XS = _airmass_for_zenith_distance(ZS)

    # Sky brightness per band (erg units internally)
    B_erg: Dict[str, float] = {}

    for i, band in enumerate(_BANDS):
        k = K[band]

        # Dark night sky baseline with 11-year solar cycle modulation
        BN = BO[i] * (1.0 + 0.3 * math.cos(6.283 * (year - 1992) / 11.0))
        BN *= (0.4 + 0.6 / math.sqrt(1.0 - 0.96 * math.sin(ZZ) ** 2))
        BN *= 10.0 ** (-0.4 * k * X)

        # Moonlight
        MM = -12.73 + 0.026 * abs(phase_moon_deg) + 4e-9 * (phase_moon_deg ** 4)
        MM += CM[i]
        C3 = 10.0 ** (-0.4 * k * XM)
        FM = 6.2e7 / (RM * RM + 1e-10) + 10.0 ** (6.15 - RM / 40.0)
        FM += 10.0 ** 5.36 * (1.06 + math.cos(_deg_to_rad(RM)) ** 2)
        BM = 10.0 ** (-0.4 * (MM - MO[i] + 43.27))
        BM *= (1.0 - 10.0 ** (-0.4 * k * X))
        BM *= (FM * C3 + 440000.0 * (1.0 - C3))

        # Twilight / daylight from Sun
        HS = 90.0 - ZS
        BT = 10.0 ** (-0.4 * (MS[i] - MO[i] + 32.5 - HS - (Z / (360.0 * k))))
        BT *= (100.0 / max(RS, 1e-6)) * (1.0 - 10.0 ** (-0.4 * k * X))

        C4 = 10.0 ** (-0.4 * k * XS)
        FS = 6.2e7 / (RS * RS + 1e-10) + 10.0 ** (6.15 - RS / 40.0)
        FS += 10.0 ** 5.36 * (1.06 + math.cos(_deg_to_rad(RS)) ** 2)
        BD = 10.0 ** (-0.4 * (MS[i] - MO[i] + 43.27))
        BD *= (1.0 - 10.0 ** (-0.4 * k * X))
        BD *= (FS * C4 + 440000.0 * (1.0 - C4))

        sky_component = BD if BD < BT else BT
        total = BN + sky_component
        if ZM < 90.0:
            total += BM

        B_erg[band] = total

    # Convert V-band brightness to nanoLamberts (nL)
    # (Explicitly listed as a corrected JS line.)
    BL_nL = B_erg["V"] / 1.02e-3

    # Threshold TH uses two regimes depending on sky brightness
    if BL_nL < 1500.0:
        C1 = 10.0 ** (-9.8)
        C2 = 10.0 ** (-1.9)
    else:
        C1 = 10.0 ** (-8.350001)
        C2 = 10.0 ** (-5.9)

    TH = C1 * (1.0 + math.sqrt(C2 * BL_nL)) ** 2

    # Limiting magnitude ("Visual Extinction Magnitude")
    # MN = -16.57 - 2.5 log10(TH) - DM(V) + 5 log10(Snellen)
    limiting_mag = -16.57 - 2.5 * math.log10(TH) - DM["V"] + 5.0 * math.log10(SN)

    # Package
    return {
        "airmass": {
            "GAS": XG,
            "AEROSOL": XA,
            "OZONE": XO,
            "STAR": X,
            "MOON": XM,
            "SUN": XS,
        },
        "extinction": {
            "k_mag_per_airmass": {b: K[b] for b in _BANDS},
            "ext_mag": {b: DM[b] for b in _BANDS},
        },
        "sky_brightness": {
            "nanolamberts_V": BL_nL,
            "picoerg": {b: B_erg[b] * 1e12 for b in _BANDS},
        },
        "threshold": TH,
        "visual_extinction_magnitude": limiting_mag,
        "geometry": {
            "TopoACRV_deg": alt_obj_deg - alt_sun_deg,
            "Sun_elongation_deg": RS,
            "Moon_elongation_deg": RM,
            "Object_altitude_deg": alt_obj_deg,
        },
    }


# -----------------------------
# Public API (web-style output)
# -----------------------------

def vislimit(
    *,
    cond: ObservingConditions,
    month: int,
    year: int,
    alt_moon_deg: float,
    az_moon_compared_to_object_deg: float,
    alt_sun_deg: float,
    az_sun_compared_to_object_deg: float,
    object_altitude_deg: float,
    goal_magnitude: Optional[float] = None,
    phase_moon_deg: Optional[float] = None,
    obs_datetime_utc: Optional[datetime] = None,
    step_deg: float = 0.01,
    max_iter: int = 20000,
) -> Dict[str, Any]:
    """Run VISLIMIT in the same two-mode style as the Reijs web page.

    Returns a dict with sections/labels that match the web calculator.
    """

    # Resolve phase
    if phase_moon_deg is None:
        if obs_datetime_utc is not None:
            phase_moon_deg = _moon_phase_angle_deg_from_datetime(obs_datetime_utc)
        else:
            # If you don't provide time, default to new moon
            phase_moon_deg = 0.0

    limiting_mode = (goal_magnitude is None) or (goal_magnitude == 99)

    resolved_inputs = {
        "Month Number": month,
        "Year": year,
        "Phase angle (deg)": phase_moon_deg,
        "Moon topo altitude (deg)": alt_moon_deg,
        "Moon azimuth compared to Object (deg)": az_moon_compared_to_object_deg,
        "Sun topo altitude (deg)": alt_sun_deg,
        "Sun azimuth compared to Object (deg)": az_sun_compared_to_object_deg,
        "Relative Humidity (%)": cond.humidity_percent,
        "Air Temperature (C)": cond.temperature_c,
        "Latitude (deg)": cond.latitude_deg,
        "Height (m)": cond.height_m,
        "Snellen Ratio": cond.snellen_ratio,
        "Object topo altitude (deg)": object_altitude_deg,
        "Goal Magnitude": goal_magnitude if goal_magnitude is not None else 99,
        "obs_datetime_utc": obs_datetime_utc.isoformat() if obs_datetime_utc is not None else None,
    }

    def eval_at(altitude_deg: float) -> Dict[str, Any]:
        return _compute_vislimit_once(
            cond=cond,
            month=month,
            year=year,
            phase_moon_deg=float(phase_moon_deg),
            alt_obj_deg=float(altitude_deg),
            alt_moon_deg=float(alt_moon_deg),
            az_moon_rel_deg=float(az_moon_compared_to_object_deg),
            alt_sun_deg=float(alt_sun_deg),
            az_sun_rel_deg=float(az_sun_compared_to_object_deg),
        )

    if limiting_mode:
        core = eval_at(object_altitude_deg)
        return {
            "inputs": resolved_inputs,
            "RESULTS": {
                "Sky Brightness: nanoLamberts": core["sky_brightness"]["nanolamberts_V"],
                "Visual Extinction Magnitude: magnitudes": core["visual_extinction_magnitude"],
                "TopoACRV: degrees": core["geometry"]["TopoACRV_deg"],
                "Sun's elongation: degrees": core["geometry"]["Sun_elongation_deg"],
            },
            "Supplemental Results": {
                "AIRMASS": core["airmass"],
                "Astromincal extinction magnitude/air mass": core["extinction"]["k_mag_per_airmass"],
                "Extinction magnitude": core["extinction"]["ext_mag"],
                "Sky Brightness picoergs/cm^2/micron/arcsec": core["sky_brightness"]["picoerg"],
                "Geometry": core["geometry"],
                "TH (threshold)": core["threshold"],
            },
            "meta": {"mode": "limiting_magnitude", "iterations": 1, "converged": True},
        }

    # altitude-search mode
    current_alt = float(object_altitude_deg)
    last: Optional[Dict[str, Any]] = None
    for it in range(1, max_iter + 1):
        last = eval_at(current_alt)
        if last["visual_extinction_magnitude"] >= float(goal_magnitude) - 1e-6:
            resolved_inputs["Object topo altitude (deg)"] = current_alt
            return {
                "inputs": resolved_inputs,
                "RESULTS": {
                    "Sky Brightness: nanoLamberts": last["sky_brightness"]["nanolamberts_V"],
                    "Visual Extinction Magnitude: magnitudes": last["visual_extinction_magnitude"],
                    "TopoACRV: degrees": last["geometry"]["TopoACRV_deg"],
                    "Sun's elongation: degrees": last["geometry"]["Sun_elongation_deg"],
                },
                "Supplemental Results": {
                    "AIRMASS": last["airmass"],
                    "Astromincal extinction magnitude/air mass": last["extinction"]["k_mag_per_airmass"],
                    "Extinction magnitude": last["extinction"]["ext_mag"],
                    "Sky Brightness picoergs/cm^2/micron/arcsec": last["sky_brightness"]["picoerg"],
                    "Geometry": last["geometry"],
                    "TH (threshold)": last["threshold"],
                },
                "meta": {"mode": "extinction_angle", "iterations": it, "converged": True},
            }
        current_alt += float(step_deg)

    # not converged
    core = last if last is not None else eval_at(current_alt)
    resolved_inputs["Object topo altitude (deg)"] = current_alt
    return {
        "inputs": resolved_inputs,
        "RESULTS": {
            "Sky Brightness: nanoLamberts": core["sky_brightness"]["nanolamberts_V"],
            "Visual Extinction Magnitude: magnitudes": core["visual_extinction_magnitude"],
            "TopoACRV: degrees": core["geometry"]["TopoACRV_deg"],
            "Sun's elongation: degrees": core["geometry"]["Sun_elongation_deg"],
        },
        "Supplemental Results": {
            "AIRMASS": core["airmass"],
            "Astromincal extinction magnitude/air mass": core["extinction"]["k_mag_per_airmass"],
            "Extinction magnitude": core["extinction"]["ext_mag"],
            "Sky Brightness picoergs/cm^2/micron/arcsec": core["sky_brightness"]["picoerg"],
            "Geometry": core["geometry"],
            "TH (threshold)": core["threshold"],
        },
        "meta": {"mode": "extinction_angle", "iterations": max_iter, "converged": False},
    }


if __name__ == "__main__":
    # Minimal smoke test, matching the two modes mentioned on the web page.
    cond = ObservingConditions(
        humidity_percent=20.0,
        temperature_c=21.0,
        latitude_deg=35.5,
        height_m=930.0,
        snellen_ratio=1.0,
    )

    out1 = vislimit(
        cond=cond,
        month=1,
        year=2024,
        alt_moon_deg=-90,
        az_moon_compared_to_object_deg=180,
        alt_sun_deg=-0.3,
        az_sun_compared_to_object_deg=0.142,
        object_altitude_deg=11.83,
        goal_magnitude=99,
        phase_moon_deg=0.0,
    )
    print("RESULTS (Goal=99)")
    for k, v in out1["RESULTS"].items():
        print(f"- {k}: {v}")

    out2 = vislimit(
        cond=cond,
        month=1,
        year=2024,
        alt_moon_deg=-90,
        az_moon_compared_to_object_deg=180,
        alt_sun_deg=-0.3,
        az_sun_compared_to_object_deg=0.142,
        object_altitude_deg=11.83,
        goal_magnitude=4.0,
        phase_moon_deg=0.0,
    )
    print("\nRESULTS (Goal=4)")
    for k, v in out2["RESULTS"].items():
        print(f"- {k}: {v}")
