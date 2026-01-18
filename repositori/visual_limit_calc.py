"""
This module provides a Python implementation of the visual limiting magnitude
calculator that appears on the Astronomical Visual Limits calculation page
published by Victor Reijs (based on work by Bradley E. Schaefer and Larry
Bogan).  The original implementation was written in JavaScript and is
available at the Archaeocosmology website.  The code here reproduces the
logic of that script in a clear and well‑structured Python form.

The core routine accepts atmospheric and observational parameters and
computes either the limiting magnitude of an object at a given altitude
or, if a goal magnitude is supplied, the altitude at which an object of
that brightness becomes just visible.  Intermediate quantities such as
air mass, extinction coefficients, sky brightness and the final visual
threshold are also returned for inspection.

References:

* The JavaScript source of the calculator defines the photometric
  constants and algorithms used to estimate extinction and sky
  brightness.  Arrays of wavelengths, zero points, ozone factors,
  aerosol factors, night sky background and colour terms are taken
  directly from that source【771736117259114†L13-L24】.  The code
  computes air masses, extinction coefficients and differential
  extinctions per waveband【771736117259114†L61-L81】, then models dark sky,
  moonlight, twilight and daylight contributions to the overall sky
  brightness【771736117259114†L90-L115】.  Finally it converts the sky
  brightness into nanolamberts and applies Schaefer’s threshold law to
  derive the visual magnitude threshold【771736117259114†L124-L135】.

Usage example:

>>> from visual_limit_calc import visual_limit
>>> result = visual_limit(
...     month=3, year=1972, phase_angle=180,  # New moon
...     altmoon=-90, azimoon=0,               # Moon below horizon
...     altsun=0, azisun=0,                   # Sun on horizon
...     humidity=20, temperature=21,          # Relative humidity and temperature
...     latitude=35.5, altitude=930,          # Observer latitude (deg) and altitude (m)
...     snellen_ratio=1.0, altstar=19.05,     # Visual acuity and object altitude
...     goal_magnitude=99                     # Compute limiting magnitude (99 indicates no goal)
... )
>>> print(result["limiting_magnitude"])
"""

from __future__ import annotations
import math
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class ExtinctionResult:
    """Container for intermediate extinction and sky brightness results."""
    sky_brightness_nl: float  # sky brightness in nanolamberts (BL)
    limiting_magnitude: float  # resulting visual limiting magnitude (MN)
    object_altitude: float  # topocentric altitude at which MN was computed
    gas_airmass: float
    aerosol_airmass: float
    ozone_airmass: float
    K: List[float] = field(default_factory=list)  # extinction coefficients per waveband
    DM: List[float] = field(default_factory=list)  # extinction/airmass per waveband
    B: List[float] = field(default_factory=list)   # sky brightness per waveband (picoergs)


def _compute_for_altitude(
    month: float,
    year: float,
    phase_angle: float,
    altmoon: float,
    azimoon: float,
    altsun: float,
    azisun: float,
    humidity: float,
    temperature: float,
    latitude: float,
    altitude: float,
    snellen_ratio: float,
    altitude_star: float,
) -> ExtinctionResult:
    """
    Compute limiting magnitude and auxiliary quantities for a fixed star altitude.

    This function implements the core of the JavaScript algorithm.  It
    calculates air mass factors, extinction coefficients (gas, aerosol,
    ozone, water vapour), sky brightness contributions (dark sky,
    moonlight, twilight and daylight) and the resultant visual threshold.
    Parameters mirror those in the web form; all angular quantities are
    expressed in degrees.

    Returns
    -------
    ExtinctionResult
        Dataclass containing the limiting magnitude, sky brightness in
        nanolamberts, input altitude, air masses and per band arrays.
    """
    # Photometric constants defined by Schaefer/Bogan【771736117259114†L13-L24】
    WA = [0.365, 0.44, 0.55, 0.7, 0.9]  # wavelengths for U,B,V,R,I bands (microns)
    MO = [-10.93, -10.45, -11.05, -11.90, -12.70]  # lunar magnitude offsets per band
    OZ = [0.000, 0.000, 0.031, 0.008, 0.000]  # ozone extinction factors
    WT = [0.074, 0.045, 0.031, 0.020, 0.015]  # water vapour coefficients
    BO = [8.0e-14, 7.0e-14, 1.0e-13, 1.0e-13, 3.0e-13]  # dark sky background (ergs)
    CM = [1.36, 0.91, 0.00, -0.76, -1.17]  # colour term corrections
    MS = [-25.96, -26.09, -26.74, -27.26, -27.55]  # solar magnitudes per band

    RD = math.pi / 180.0  # degrees to radians conversion

    # Convert azimuth/altitude differences to elongations (angular separations)
    # RM: angle between Moon and object; RS: angle between Sun and object【771736117259114†L51-L55】
    RM = math.degrees(
        math.acos(
            math.sin(altmoon * RD) * math.sin(altitude_star * RD)
            + math.cos(altmoon * RD) * math.cos(altitude_star * RD) * math.cos(azimoon * RD)
        )
    )
    RS = math.degrees(
        math.acos(
            math.sin(altsun * RD) * math.sin(altitude_star * RD)
            + math.cos(altsun * RD) * math.cos(altitude_star * RD) * math.cos(azisun * RD)
        )
    )

    # Zenith distances【771736117259114†L56-L60】
    ZM = 90.0 - altmoon
    ZS = 90.0 - altsun
    Z = 90.0 - altitude_star

    # Extinction subroutine【771736117259114†L61-L81】
    LT = latitude * RD  # latitude in radians
    RA = (month - 3.0) * 30.0 * RD  # right ascension proxy (radians)
    # Determine sign of latitude; avoid division by zero
    SL = 1.0 if latitude >= 0 else -1.0
    ZZ = Z * RD
    # Airmass calculations【771736117259114†L61-L81】
    XG = 1.0 / (math.cos(ZZ) + 0.0286 * math.exp(-10.5 * math.cos(ZZ)))
    XA = 1.0 / (math.cos(ZZ) + 0.0123 * math.exp(-24.5 * math.cos(ZZ)))
    # Ozone airmass uses Earth radius and 20 km ozone layer height【771736117259114†L66-L70】
    sinZZ = math.sin(ZZ)
    XO = 1.0 / math.sqrt(1.0 - (sinZZ / (1.0 + (20.0 / 6378.0))) ** 2)

    K: List[float] = [0.0] * 5  # extinction coefficients per band
    DM: List[float] = [0.0] * 5  # differential magnitudes per band
    B: List[float] = [0.0] * 5   # sky brightness per band (picoergs)

    # Loop over photometric bands to compute extinction coefficients and DM【771736117259114†L72-L81】
    for i in range(5):
        KR = 0.1066 * math.exp(-altitude / 8200.0) * (WA[i] / 0.55) ** (-4.0)
        KA = 0.1 * (WA[i] / 0.55) ** (-1.3) * math.exp(-altitude / 1500.0)
        # Humidity and seasonal adjustment with corrected 0.33 factor【771736117259114†L72-L80】
        # Protect against log(0) by enforcing minimum humidity (>= 1%).
        rh_fraction = max(humidity / 100.0, 0.01)
        KA *= (1.0 - 0.32 / math.log(rh_fraction)) ** 1.33
        KA *= (1.0 + 0.33 * SL * math.sin(RA))
        KO = OZ[i] * (3.0 + 0.4 * (LT * math.cos(RA) - math.cos(3.0 * LT))) / 3.0
        KW = WT[i] * 0.94 * (humidity / 100.0) * math.exp(temperature / 15.0) * math.exp(-altitude / 8200.0)
        K[i] = KR + KA + KO + KW
        DM[i] = KR * XG + KA * XA + KO * XO + KW * XG

    # Sky brightness subroutine【771736117259114†L90-L120】
    X = 1.0 / (math.cos(ZZ) + 0.025 * math.exp(-11.0 * math.cos(ZZ)))
    # Air mass of Moon and Sun for scattering calculations
    # Zenith distance of Moon and Sun in radians
    XM = 1.0 / (math.cos(ZM * RD) + 0.025 * math.exp(-11.0 * math.cos(ZM * RD)))
    if ZM > 90.0:
        XM = 40.0
    XS = 1.0 / (math.cos(ZS * RD) + 0.025 * math.exp(-11.0 * math.cos(ZS * RD)))
    if ZS > 90.0:
        XS = 40.0

    for i in range(5):
        # Dark night sky brightness【771736117259114†L90-L93】
        BN = BO[i] * (1.0 + 0.3 * math.cos(6.283 * (year - 1992.0) / 11.0))
        BN *= (0.4 + 0.6 / math.sqrt(max(1.0 - 0.96 * (math.sin(ZZ) ** 2), 1e-6)))
        BN *= 10.0 ** (-0.4 * K[i] * X)
        # Moonlight brightness【771736117259114†L94-L103】
        MM = -12.73 + 0.026 * abs(phase_angle) + 4e-9 * (phase_angle ** 4)
        MM += CM[i]
        C3 = 10.0 ** (-0.4 * K[i] * XM)
        # Avoid division by zero in RM when Moon is directly overhead
        if RM == 0.0:
            FM = 0.0
        else:
            FM = (6.2e7) / (RM ** 2) + 10.0 ** (6.15 - RM / 40.0)
        FM += 10.0 ** 5.36 * (1.06 + (math.cos(RM * RD) ** 2))
        BM = 10.0 ** (-0.4 * (MM - MO[i] + 43.27))
        BM *= (1.0 - 10.0 ** (-0.4 * K[i] * X))
        BM *= (FM * C3 + 440000.0 * (1.0 - C3))
        # Twilight brightness【771736117259114†L104-L107】
        HS = 90.0 - ZS
        # Avoid division by zero in K[i] and RS
        if K[i] == 0.0 or RS == 0.0:
            BT = 0.0
        else:
            BT = 10.0 ** (-0.4 * (MS[i] - MO[i] + 32.5 - HS - (Z / (360.0 * K[i]))))
            BT *= (100.0 / RS) * (1.0 - 10.0 ** (-0.4 * K[i] * X))
        # Daylight brightness【771736117259114†L108-L115】
        C4 = 10.0 ** (-0.4 * K[i] * XS)
        # Avoid division by zero in RS again
        if RS == 0.0:
            FS = 0.0
        else:
            FS = 6.2e7 / (RS ** 2) + 10.0 ** (6.15 - RS / 40.0)
        FS += 10.0 ** 5.36 * (1.06 + (math.cos(RS * RD) ** 2))
        BD = 10.0 ** (-0.4 * (MS[i] - MO[i] + 43.27))
        BD *= (1.0 - 10.0 ** (-0.4 * K[i] * X))
        BD *= (FS * C4 + 440000.0 * (1.0 - C4))
        # Combine sky brightness components【771736117259114†L115-L120】
        # Use the smaller of daylight and twilight contributions
        if BD < BT:
            brightness = BN + BD
        else:
            brightness = BN + BT
        # Add moonlight contribution if Moon is above horizon
        if ZM < 90.0:
            brightness += BM
        # Convert from ergs to picoergs
        B[i] = brightness * 1.0e12

    # Convert V band (index 2) sky brightness to nanolamberts【771736117259114†L124-L135】
    # The factor 1.02e-3 converts picoergs/cm²/μm/arcsec² to nanolamberts (approximate)
    BL = B[2] / 1.02e-3
    # Determine empirical constants for visual threshold【771736117259114†L124-L135】
    if BL < 1500.0:
        C1 = 10.0 ** (-9.8)
        C2 = 10.0 ** (-1.9)
    else:
        C1 = 10.0 ** (-8.350001)
        C2 = 10.0 ** (-5.9)
    TH = C1 * (1.0 + math.sqrt(C2 * BL)) ** 2.0
    # Limiting magnitude incorporating extinction and eye sensitivity【771736117259114†L135-L139】
    # Use math.log10 for base‑10 logarithms
    limiting_mag = -16.57 - 2.5 * math.log10(TH) - DM[2] + 5.0 * math.log10(max(snellen_ratio, 1e-6))
    return ExtinctionResult(
        sky_brightness_nl=BL,
        limiting_magnitude=limiting_mag,
        object_altitude=altitude_star,
        gas_airmass=XG,
        aerosol_airmass=XA,
        ozone_airmass=XO,
        K=K,
        DM=DM,
        B=B,
    )


def visual_limit(
    month: float,
    year: float,
    phase_angle: float,
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
    goal_magnitude: Optional[float] = None,
    step_deg: float = 0.01,
) -> Dict[str, Any]:
    """
    Compute the visual limiting magnitude or the altitude required to reach
    a given magnitude.

    Parameters follow those of the original JavaScript form.  Angles are
    given in degrees.  If `goal_magnitude` is set to 99 or left as
    ``None``, the function computes the limiting magnitude for the
    supplied star altitude.  Otherwise it will iterate starting at
    altitude 0° and increment by `step_deg` until the computed
    magnitude equals or exceeds the goal magnitude.  A small step
    produces finer resolution at the cost of more computation.

    Returns
    -------
    Dict[str, Any]
        A dictionary containing the sky brightness (nanolamberts), the
        limiting magnitude, the altitude at which the result was
        obtained, air mass components and the per‑band extinction and
        brightness arrays.  Fields include:

        - ``sky_brightness`` – sky brightness in nanolamberts
        - ``limiting_magnitude`` – limiting magnitude of the object
        - ``altitude`` – the altitude used/found (degrees)
        - ``gas_airmass`` – air mass for gas scattering
        - ``aerosol_airmass`` – air mass for aerosols
        - ``ozone_airmass`` – air mass for ozone
        - ``K`` – list of extinction coefficients per photometric band
        - ``DM`` – list of differential magnitudes per band
        - ``B`` – list of sky brightness per band (picoergs)
    """
    # Set sentinel value according to original script
    if goal_magnitude is None:
        goal = 99.0
    else:
        goal = float(goal_magnitude)

    # Starting altitude: if computing limiting magnitude, start at the provided star altitude
    # Otherwise start searching from horizon (0°).
    if goal == 99.0:
        start_alt = altstar
    else:
        start_alt = 0.0

    MN = -99.0  # initialise limiting magnitude to a very faint value
    current_alt = start_alt
    last_result: Optional[ExtinctionResult] = None

    # Iterate until the limiting magnitude reaches/exceeds the goal
    while MN < goal:
        result = _compute_for_altitude(
            month,
            year,
            phase_angle,
            altmoon,
            azimoon,
            altsun,
            azisun,
            humidity,
            temperature,
            latitude,
            altitude,
            snellen_ratio,
            current_alt,
        )
        MN = result.limiting_magnitude
        last_result = result
        # If computing for a fixed altitude, copy the computed magnitude into the goal
        # so that the loop terminates after a single iteration (mimicking the
        # JavaScript behaviour where GoalMg==99 triggers GoalMg = MN in the loop)【771736117259114†L140-L144】.
        if goal == 99.0:
            goal = MN
        else:
            # Continue searching: increment altitude by specified step if threshold not yet reached
            if MN < goal:
                current_alt += step_deg
    # Format the result dictionary for user consumption
    if last_result is None:
        # Should never happen; return empty result in pathological case
        return {}
    return {
        "sky_brightness": last_result.sky_brightness_nl,
        "limiting_magnitude": last_result.limiting_magnitude,
        "altitude": last_result.object_altitude,
        "gas_airmass": last_result.gas_airmass,
        "aerosol_airmass": last_result.aerosol_airmass,
        "ozone_airmass": last_result.ozone_airmass,
        "K": last_result.K,
        "DM": last_result.DM,
        "B": last_result.B,
    }


__all__ = ["visual_limit"]