"""
Implementation of the correction factors used in Bradley E. Schaefer's
"Telescopic Limiting Magnitudes" model.

This module collects each of the numbered equations from the paper and
exposes them as Python functions.  The aim is to allow a caller to
compute the limiting magnitude of a telescope–eye system under a given
set of observing conditions by following Schaefer's algorithm.  All of
the expressions implemented here come from the paper "Telescopic
Limiting Magnitudes" (see the uploaded PDF for full context).

The core quantities are:

* ``I``: the basic eye sensitivity to a point source against a uniform
  background (Equation 2).  This depends on the sky brightness ``B``
  expressed in millimicro‑Lamberts (mμL) and the appropriate
  constants ``C`` and ``K`` for photopic or scotopic vision.

* ``F_e``: the atmospheric extinction factor (Equation 3).

* ``F_t``: the transmission factor of the telescope including losses
  through optical surfaces and any central obstruction (Equation 4).

* ``D_e``: the observer's eye pupil diameter as a function of age
  (Equation 5).

* ``F_p``: the pupil factor comparing the telescope's exit pupil to
  the observer's eye pupil (Equation 6).

* ``F_a``: the area scaling factor when the eye pupil is smaller than
  the telescope aperture (additional unnumbered expression in the
  paper).

* ``F_m``: the magnification factor, which accounts for the increase
  in apparent sky brightness due to magnification (additional
  unnumbered expression in the paper).

* ``F_r``: the seeing factor that accounts for seeing disc blurring
  (Equation 7).

* ``F_sc``: the Stiles–Crawford efficiency correction (Equation 9).

* ``F_c``: the colour correction for the difference between the
  photopic and scotopic passbands and the star's (B−V) colour index
  (Equation 13).

After computing these factors one can combine them as in
Equation (14) to obtain the effective irradiance ``I*`` on the retina,
and finally convert to a limiting magnitude ``m`` via Equation (16).

This file exposes helper functions to compute each term and a
convenience routine ``compute_limiting_magnitude`` that returns the
final magnitude directly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def eye_sensitivity(B: float) -> float:
    """
    Compute the base eye sensitivity ``I`` as a function of sky
    brightness ``B`` using Equation (2).

    Parameters
    ----------
    B : float
        Sky background brightness in millimicro‑Lamberts (mμL).

    Returns
    -------
    float
        The brightness ``I`` of a just‑detectable point source in
        foot‑candles.  The constants ``C`` and ``K`` are chosen
        according to whether the regime is photopic (log(B) > 3.17) or
        scotopic (log(B) < 3.17).
    """
    if B <= 0:
        raise ValueError("Sky brightness B must be positive")
    logB = math.log10(B)
    if logB < 3.17:
        # Scotopic regime
        C = 10 ** (-9.80)
        K = 10 ** (-1.90)
    else:
        # Photopic regime
        C = 10 ** (-8.35)
        K = 10 ** (-5.90)
    return C * (1.0 + math.sqrt(K * B)) ** 2


def atmospheric_extinction_factor(B: float, k_V: float, Z_deg: float) -> float:
    """
    Compute the atmospheric extinction correction factor ``F_e``
    according to Equation (3).

    The paper expresses the correction through a logarithmic relation
    ``2.5 log10(F_e) = q * k_V * sec(Z)``.  We solve this for ``F_e``.

    Parameters
    ----------
    B : float
        Sky brightness in mμL (used only to choose the ``q`` constant).
    k_V : float
        V‑band extinction coefficient (magnitudes per air mass).  A
        typical value in clear skies is ~0.2–0.3 mag.
    Z_deg : float
        Zenith distance in degrees.  ``Z = 0`` corresponds to the
        zenith and ``Z = 90`` to the horizon.

    Returns
    -------
    float
        ``F_e``, the factor by which the limiting irradiance increases
        due to atmospheric extinction.
    """
    if B <= 0:
        raise ValueError("Sky brightness B must be positive")
    # Choose q depending on whether the eye operates in the scotopic
    # (log B < 3.17) or photopic (log B > 3.17) regime.
    logB = math.log10(B)
    q = 1.2 if logB < 3.17 else 1.0
    # Convert zenith distance to air mass using secant of Z.  We
    # protect against cos(Z) = 0 near the horizon by clamping Z.
    Z_rad = math.radians(Z_deg)
    # Avoid division by zero; restrict Z to < 89.5° to keep sec(Z) finite.
    cosZ = math.cos(Z_rad)
    if cosZ <= 0:
        # At the horizon the air mass formally diverges; use a large
        # number to approximate extremely high extinction.  This is
        # physically unlikely for visible observations.
        secZ = 1e6
    else:
        secZ = 1.0 / cosZ
    # Solve 2.5 log10(F_e) = q * k_V * sec(Z)
    exponent = q * k_V * secZ / 2.5
    return 10 ** exponent


def transmission_factor(t: float, n: int, D_s: float, D: float) -> float:
    """
    Compute the telescope transmission factor ``F_t`` from
    Equation (4).

    ``F_t = 1 / [t**n * (1 − (D_s/D)**2)]``

    Parameters
    ----------
    t : float
        Average transmission of a single optical surface (0 < t ≤ 1).
    n : int
        Number of optical surfaces through which light passes.  For
        example, a refractor with an objective and eyepiece has
        ``n=4``.
    D_s : float
        Diameter of the central obstruction (same units as ``D``).
    D : float
        Telescope aperture diameter (same units as ``D_s``).

    Returns
    -------
    float
        ``F_t`` – the multiplicative factor by which the irradiance
        increases to compensate for optical losses.
    """
    if not (0 < t <= 1):
        raise ValueError("t must be between 0 and 1")
    if D <= 0:
        raise ValueError("D (aperture) must be positive")
    if D_s < 0:
        raise ValueError("D_s (obstruction diameter) must be non‑negative")
    if n < 0:
        raise ValueError("n (number of surfaces) must be non‑negative")
    # Transmission through n surfaces
    trans = t ** n
    # Central obstruction area factor
    if D_s >= D:
        obstruction_factor = 0.0
    else:
        obstruction_factor = 1.0 - (D_s / D) ** 2
    if obstruction_factor <= 0:
        raise ValueError(
            "Central obstruction diameter must be smaller than aperture to obtain a non‑zero transmission factor"
        )
    return 1.0 / (trans * obstruction_factor)


def eye_pupil_diameter(age: float) -> float:
    """
    Compute the observer's pupil diameter ``D_e`` in millimetres
    according to Equation (5).

    ``D_e = 7 mm * exp( −0.5 * (age / 100)**2 )``

    Parameters
    ----------
    age : float
        Observer's age in years.

    Returns
    -------
    float
        Pupil diameter in millimetres.
    """
    if age < 0:
        raise ValueError("age must be non‑negative")
    exponent = -0.5 * (age / 100.0) ** 2
    return 7.0 * math.exp(exponent)


def exit_pupil(D: float, M: float) -> float:
    """
    Compute the telescope's exit pupil diameter ``D_exit``.

    The exit pupil is simply the aperture divided by the magnification
    (``D_exit = D / M``).  ``D`` and ``M`` must be in the same
    units for ``D_exit`` to have the same units as ``D``.

    Parameters
    ----------
    D : float
        Telescope aperture diameter.
    M : float
        Telescope magnification.

    Returns
    -------
    float
        Exit pupil diameter.
    """
    if D <= 0:
        raise ValueError("D must be positive")
    if M <= 0:
        raise ValueError("M must be positive")
    return D / M


def pupil_factor(D: float, M: float, age: float) -> float:
    """
    Compute the pupil factor ``F_p`` from Equation (6).

    ``F_p = (D / (M * D_e))**2`` if the exit pupil (``D/M``) is
    larger than the observer's pupil (``D_e``); otherwise ``F_p = 1``.

    Parameters
    ----------
    D : float
        Telescope aperture diameter.
    M : float
        Telescope magnification.
    age : float
        Observer's age in years.

    Returns
    -------
    float
        ``F_p``, the factor needed when the eye pupil is smaller than
        the exit pupil.
    """
    D_e = eye_pupil_diameter(age)
    D_exit = exit_pupil(D, M)
    if D_e < D_exit:
        return (D / (M * D_e)) ** 2
    else:
        return 1.0


def area_factor(D: float, age: float) -> float:
    """
    Compute the area scaling factor ``F_a``.

    According to the unnumbered expressions following Equation (6) in
    the paper, ``F_a = (D_e / D)**2``.  This factor applies when the
    eye pupil ``D_e`` is smaller than the telescope aperture ``D``.

    Parameters
    ----------
    D : float
        Telescope aperture diameter.
    age : float
        Observer's age.

    Returns
    -------
    float
        ``F_a``, the area correction factor due to the eye pupil.
    """
    D_e = eye_pupil_diameter(age)
    if D <= 0:
        raise ValueError("D must be positive")
    return (D_e / D) ** 2


def magnification_factor(M: float) -> float:
    """
    Compute the magnification factor ``F_m``.

    In the unnumbered expressions following Equation (6) the paper
    defines ``F_m = M**2``.  This factor represents the increase in
    apparent sky brightness in the eyepiece as magnification grows.

    Parameters
    ----------
    M : float
        Telescope magnification.

    Returns
    -------
    float
        ``F_m = M**2``.
    """
    if M <= 0:
        raise ValueError("Magnification M must be positive")
    return M ** 2


def seeing_factor(theta_arcsec: float, M: float) -> float:
    """
    Compute the seeing correction factor ``F_r`` according to
    Equation (7).

    ``F_r = sqrt( (2 * theta * M) / 900 )`` if ``2 * theta * M > 900``;
    otherwise ``F_r = 1``.

    ``theta`` is the seeing disk radius in arcseconds, and ``M`` is
    magnification.  The constant ``900`` comes from the paper.

    Parameters
    ----------
    theta_arcsec : float
        Radius of the seeing disc in arcseconds.
    M : float
        Telescope magnification.

    Returns
    -------
    float
        ``F_r``, the seeing correction factor.
    """
    if theta_arcsec < 0:
        raise ValueError("theta_arcsec must be non‑negative")
    if M <= 0:
        raise ValueError("M must be positive")
    # Determine whether the star image becomes an extended source.
    # 2 * theta * M is the apparent diameter; if it exceeds 900 arcsec
    # (15 arcmin) the star is extended and the detection threshold
    # improves proportionally to the square root of the area ratio.
    product = 2.0 * theta_arcsec * M
    if product > 900.0:
        return math.sqrt(product / 900.0)
    else:
        return 1.0


def stiles_crawford_factor(B: float, D: float, M: float, age: float) -> float:
    """
    Compute the Stiles–Crawford correction factor ``F_sc`` from
    Equation (9).

    The Stiles–Crawford effect reduces the efficiency of light rays
    entering the pupil away from its centre.  The correction depends on
    whether the regime is photopic (log10(B) > 3.17) or scotopic
    (log10(B) < 3.17) and on the ratio of the telescope's exit pupil
    to the observer's eye pupil.

    When ``log10(B) > 3.17`` (photopic), the factor is:
    ``F_sc = (D_e * M / D) * [1 − exp(−0.026 * (D/M)**2)] / [1 − exp(−0.026 * D_e**2)]``

    When ``log10(B) < 3.17`` (scotopic), the factor is:
    ``F_sc = [1 − (D / (12.4 * M))**4] / [1 − (D_e / 12.4)**4]``

    These expressions apply only when the eye pupil is larger than the
    exit pupil; otherwise the factor is taken as 1.

    Parameters
    ----------
    B : float
        Sky brightness in mμL.
    D : float
        Telescope aperture diameter.
    M : float
        Telescope magnification.
    age : float
        Observer's age.

    Returns
    -------
    float
        ``F_sc``.  Returns ``1`` if the observer's eye pupil is
        smaller than or equal to the exit pupil, in which case rays do
        not sample significantly different parts of the pupil.
    """
    if B <= 0:
        raise ValueError("Sky brightness B must be positive")
    if D <= 0:
        raise ValueError("Aperture D must be positive")
    if M <= 0:
        raise ValueError("Magnification M must be positive")

    D_e = eye_pupil_diameter(age)
    exit_pup = exit_pupil(D, M)
    # The correction is applied only when the eye pupil exceeds the exit pupil.
    if D_e <= exit_pup:
        return 1.0

    logB = math.log10(B)
    if logB > 3.17:
        # Photopic regime
        numerator = D_e * M / D
        # Avoid overflow/underflow in exponentials for large arguments.
        exp_term_exit = math.exp(-0.026 * (D / M) ** 2)
        exp_term_eye = math.exp(-0.026 * D_e ** 2)
        # Denominator cannot be zero because exp_term_eye is always < 1.
        F_sc = numerator * (1.0 - exp_term_exit) / (1.0 - exp_term_eye)
    else:
        # Scotopic regime
        base_exit = D / (12.4 * M)
        base_eye = D_e / 12.4
        # Avoid negative values if either term exceeds unity.
        one_minus_exit4 = 1.0 - min(base_exit ** 4, 1.0)
        one_minus_eye4 = 1.0 - min(base_eye ** 4, 1.0)
        # If the denominator is zero (both terms saturate), the factor
        # tends to 1 because both numerator and denominator go to zero.
        if abs(one_minus_eye4) < 1e-12:
            F_sc = 1.0
        else:
            F_sc = one_minus_exit4 / one_minus_eye4
    # Constrain the factor to be positive.
    return max(F_sc, 0.0)


def colour_correction_factor(B: float, B_minus_V: float) -> float:
    """
    Compute the colour correction factor ``F_c`` using Equation (13).

    The paper gives ``−2.5 log10(F_c) = 1 − (B−V)/2`` in the scotopic
    regime (``log10(B) < 3.17``) and ``−2.5 log10(F_c) = 0`` in the
    photopic regime.  We solve for ``F_c``.

    Parameters
    ----------
    B : float
        Sky brightness in mμL.
    B_minus_V : float
        (B−V) colour index of the target star.  Typical values are
        around 0 for an A0 star and 1.0 for a cool red star.

    Returns
    -------
    float
        ``F_c``.  If the regime is photopic the factor is 1; if
        scotopic then ``F_c = 10**( −(1 − (B_minus_V)/2) / 2.5 )``.
    """
    if B <= 0:
        raise ValueError("Sky brightness B must be positive")
    logB = math.log10(B)
    if logB < 3.17:
        exponent = -(1.0 - (B_minus_V) / 2.0) / 2.5
        return 10 ** exponent
    else:
        return 1.0


def sky_brightness_from_mag(B_V: float) -> float:
    """
    Convert a sky surface brightness from magnitudes per square
    arcsecond to millimicro‑Lamberts (mμL).

    The conversion comes from Equation (17):

    ``B_s = 34.08 * exp( 20.7233 - 0.92104 * B_V )``

    where ``B_V`` is the sky brightness in visual magnitudes per
    square arcsecond.

    Parameters
    ----------
    B_V : float
        Sky brightness in mag/arcsec².

    Returns
    -------
    float
        Sky brightness ``B_s`` in millimicro‑Lamberts.
    """
    exponent = 20.7233 - 0.92104 * B_V
    return 34.08 * math.exp(exponent)


def limiting_magnitude(I_star: float) -> float:
    """
    Convert the effective irradiance ``I*`` to a limiting magnitude
    using Equation (16).

    ``m = −16.57 − 2.5 log10(I*)``

    Parameters
    ----------
    I_star : float
        Effective irradiance on the retina (foot‑candles).

    Returns
    -------
    float
        Limiting magnitude ``m``.
    """
    if I_star <= 0:
        raise ValueError("I_star must be positive")
    return -16.57 - 2.5 * math.log10(I_star)


@dataclass
class SchaeferInputs:
    """
    Bundle together the parameters required to compute the limiting
    magnitude using Schaefer's method.

    Attributes
    ----------
    sky_brightness_mmag : float
        Sky background brightness in mag/arcsec².  Use
        ``sky_brightness_from_mag`` to convert to ``B`` in mμL.
    aperture : float
        Telescope aperture diameter (mm or another consistent unit).
    central_obstruction : float
        Diameter of the telescope's central obstruction in the same
        units as ``aperture`` (0 if unobstructed).
    transmission_per_surface : float
        Average transmission per optical surface (0 < t ≤ 1).
    n_surfaces : int
        Number of optical surfaces (objective elements + eyepiece
        elements).  For example, a simple refractor might have two
        lenses in the objective and two in the eyepiece (``n=4``).
    magnification : float
        Telescope magnification ``M``.
    seeing : float
        Seeing disc radius ``theta`` in arcseconds.
    extinction_coefficient : float
        Atmospheric extinction coefficient ``k_V`` in magnitudes per
        air mass.
    zenith_distance : float
        Zenith distance ``Z`` in degrees (0 at zenith, 90 at horizon).
    age : float
        Observer's age in years (used to compute eye pupil diameter).
    colour_index : float
        (B−V) colour index of the star.  Typical values: 0 for blue
        stars, 0.6 for the Sun, 1.0 for red giants.
    F_b : float
        Background scaling factor.  In Schaefer's derivation this
        factor accounts for differences between measured and modelled
        sky brightness; if unsure, use ``1.0``.
    F_s : float
        Optional additional scaling factor accounting for miscellaneous
        empirical effects (e.g. observer skill, hyperventilation).  Set
        to 1.0 if unknown.  You can instead apply the empirical
        corrections separately using the functions ``experience_correction``
        and ``hyperventilation_correction`` defined below.
    """

    sky_brightness_mmag: float
    aperture: float
    central_obstruction: float
    transmission_per_surface: float
    n_surfaces: int
    magnification: float
    seeing: float
    extinction_coefficient: float
    zenith_distance: float
    age: float
    colour_index: float
    F_b: float = 1.0
    F_s: float = 1.0


def compute_limiting_magnitude(inputs: SchaeferInputs) -> float:
    """
    Compute the limiting magnitude for the given observing setup.

    This function follows Schaefer's model closely: it converts the
    sky brightness from mag/arcsec² to mμL, determines the base eye
    sensitivity ``I``, applies all of the correction factors, and
    finally converts the resulting irradiance to a magnitude.  The
    factors ``F_b`` and ``F_s`` can be provided directly via the
    ``SchaeferInputs`` dataclass; they default to 1.0.

    Parameters
    ----------
    inputs : SchaeferInputs
        The input parameters describing the telescope, the observer,
        and the observing conditions.

    Returns
    -------
    float
        The limiting magnitude of a star that can be detected under
        the specified conditions.
    """
    # Convert the sky brightness from mag/arcsec² to B in mμL.
    B = sky_brightness_from_mag(inputs.sky_brightness_mmag)
    # Base eye sensitivity (Equation 2)
    I = eye_sensitivity(B)
    # Atmospheric extinction (Equation 3)
    F_e = atmospheric_extinction_factor(B, inputs.extinction_coefficient, inputs.zenith_distance)
    # Telescope transmission and obstruction (Equation 4)
    F_t = transmission_factor(inputs.transmission_per_surface, inputs.n_surfaces,
                               inputs.central_obstruction, inputs.aperture)
    # Pupil factor (Equation 6)
    F_p = pupil_factor(inputs.aperture, inputs.magnification, inputs.age)
    # Area factor (additional expression following Equation 6)
    F_a = area_factor(inputs.aperture, inputs.age)
    # Seeing correction (Equation 7)
    F_r = seeing_factor(inputs.seeing, inputs.magnification)
    # Stiles–Crawford correction (Equation 9)
    F_sc = stiles_crawford_factor(B, inputs.aperture, inputs.magnification, inputs.age)
    # Colour correction (Equation 13)
    F_c = colour_correction_factor(B, inputs.colour_index)
    # Magnification factor for background (un‐numbered expression)
    F_m = magnification_factor(inputs.magnification)
    # Combine all factors to get I* (Equation 14)
    I_star = I * inputs.F_b * F_e * F_t * F_p * F_a * F_r * F_sc * F_c * inputs.F_s
    # Convert to magnitude (Equation 16)
    m = limiting_magnitude(I_star)
    return m






def experience_correction(m: float, experience: float) -> float:
    """
    Apply the empirical observer experience correction to a computed
    limiting magnitude.

    According to Equation (20), ``m_corrected = m + 0.16 * (experience − 6)``
    where ``experience`` is the number of years of observing practice.

    Parameters
    ----------
    m : float
        Limiting magnitude computed without experience correction.
    experience : float
        The observer's experience in years.  An average amateur
        observer might have ``experience=6``.

    Returns
    -------
    float
        The corrected limiting magnitude.
    """
    return m + 0.16 * (experience - 6.0)


def hyperventilation_correction(m: float) -> float:
    """
    Apply the hyperventilation correction to a computed limiting
    magnitude.

    Equation (21) suggests that deliberate hyperventilation can
    improve the limiting magnitude by approximately 0.3 mag.

    Parameters
    ----------
    m : float
        Limiting magnitude before applying the hyperventilation correction.

    Returns
    -------
    float
        Limiting magnitude after adding 0.3 mag.
    """
    return m + 0.3


__all__ = [
    "SchaeferInputs",
    "compute_limiting_magnitude",
    "eye_sensitivity",
    "atmospheric_extinction_factor",
    "transmission_factor",
    "eye_pupil_diameter",
    "exit_pupil",
    "pupil_factor",
    "area_factor",
    "magnification_factor",
    "seeing_factor",
    "stiles_crawford_factor",
    "colour_correction_factor",
    "sky_brightness_from_mag",
    "limiting_magnitude",
    "experience_correction",
    "hyperventilation_correction",
]