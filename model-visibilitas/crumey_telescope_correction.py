
"""
Crumey (2014) contrast-threshold + telescope correction for extended targets (like the crescent).

Core idea:
- A telescope does NOT magically multiply the object's contrast vs sky by M^2.
- For extended sources, both object luminance and sky luminance at the eye scale similarly.
- The improvement comes from:
  1) increased apparent target area (A_a = A * M^2)
  2) changed adaptation/background luminance at the eye (B_a), determined by exit pupil and throughput

This module implements:
- Crumey's general threshold contrast for finite targets (Eq. 41)
- Crumey's telescopic threshold curve (Eq. 77) using B_a (Eq. 66) and A_a

Inputs:
- L_nl, B_nl : luminance of crescent and sky in nanolamberts (nL)
- A_arcmin2  : angular area of the crescent (arcmin^2)
- Telescope spec: aperture, magnification, throughput (transmission, surfaces, obstruction), observer age

Outputs:
- C_obj (Weber contrast) and a visibility "margin": log10(|C_obj| / C_thr)
  margin > 0 => predicted detectable; margin < 0 => not detectable.

Important caveats:
- Crumey notes he used only *positive-contrast* Blackwell data; for negative contrast you can use abs(C)
  as an approximation but expect systematic error (calibrate with observations).
- Crumey also notes the simple scotopic model is appropriate up to about 0.1 cd/m^2 (~15 mag/arcsec^2),
  but the full threshold model spans wider B using q(B). See citations in your PDF.

Reference:
- Crumey, A. (2014). Human contrast threshold and astronomical visibility. (uploaded as crumey.pdf)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Tuple

from kastner import crescent_area


# ------------------- Crescent area helper -------------------
def crescent_area_arcmin2(elongation_deg: float, r_deg: float) -> float:
    """
    Menghitung luas sabit bulan dalam arcmin² menggunakan fungsi crescent_area dari kastner.

    Parameter:
    elongation_deg (float): Sudut elongasi (derajat)
    r_deg (float): Semidiameter bulan (derajat)

    Returns:
    float: Luas sabit bulan (arcmin²)

    Catatan:
    crescent_area() mengembalikan nilai dalam derajat², dikonversi ke arcmin² (x 3600)
    """
    A_deg2 = crescent_area(elongation_deg, r_deg)
    return A_deg2 * 3600.0  # 1 derajat² = 3600 arcmin²


# Unit conversion between Schaefer-style nanolambert and Crumey's cd/m^2
# 1 lambert = 10000/pi cd/m^2 ; 1 nL = 1e-9 lambert
NL_TO_CD_M2 = (10000.0 / math.pi) * 1e-9
CD_M2_TO_NL = 1.0 / NL_TO_CD_M2


@dataclass(frozen=True)
class TelescopeSpec:
    aperture_mm: float
    magnification: float
    obstruction_mm: float = 0.0               # central obstruction diameter (mm), 0 for refractor
    transmission_per_surface: float = 0.98    # per-surface transmission/reflection
    n_surfaces: int = 4                       # number of surfaces (mirrors+glass etc.)
    extra_throughput: float = 1.0             # optional extra factor (filters, etc.)
    observer_age: float = 22.0

    # Conservative multipliers:
    # - field_factor: raises threshold to reflect non-ideal shapes/seeing; Crumey mentions "common-sense"
    #   visibility vs lab data can differ by ~2.4 (example in his paper).
    field_factor: float = 1.5

    # - phi: additional multiplier used in Crumey's telescopic threshold curve (Eq. 77).
    #   If you don't know what to set, keep 1.0 (conservative).
    phi: float = 1.0


def eye_pupil_diameter_mm(age: float) -> float:
    """
    Same functional form as in your limit_teleskop.py.
    Returns pupil diameter in mm.
    """
    age = max(0.0, min(120.0, float(age)))
    return 7.0 - (5.0 / 3.0) * (1.0 - math.exp(-0.5 * (age / 100.0) ** 2))


def telescope_throughput(spec: TelescopeSpec) -> float:
    """
    Total physical throughput T (0..1):
      T = t^n * (1 - (d/D)^2) * extra_throughput
    """
    D = float(spec.aperture_mm)
    d = float(spec.obstruction_mm)
    if D <= 0:
        raise ValueError("aperture_mm must be > 0")
    if d < 0 or d >= D:
        raise ValueError("obstruction_mm must satisfy 0 <= d < aperture_mm")

    t = float(spec.transmission_per_surface)
    n = int(spec.n_surfaces)
    if not (0 < t <= 1):
        raise ValueError("transmission_per_surface must be in (0,1]")
    if n < 0:
        raise ValueError("n_surfaces must be >= 0")

    T = (t ** n) * (1.0 - (d / D) ** 2) * float(spec.extra_throughput)
    return max(0.0, min(1.0, T))


# ------------------- Crumey threshold model -------------------
# q(B): Eq. 42-44 in your PDF:
# q = 0.6 (B < 0.193)
# q = 0.8861 + 0.4 log10(B) (0.193 <= B <= 3.426)
# q = 1.146 - 0.0885 log10(B) (B > 3.426)
_B_SCOT = 0.193
_B_PHOT = 3.426

def q_of_B(B_cd_m2: float) -> float:
    B = float(B_cd_m2)
    if B <= 0:
        return 0.6
    if B < _B_SCOT:
        return 0.6
    lg = math.log10(B)
    if B <= _B_PHOT:
        return 0.8861 + 0.4 * lg
    return 1.146 - 0.0885 * lg


# Full-curve R(B): Eq. 28 with coefficients from your PDF:
# R = ( sqrt(a1 B^-1/2 + a2 B^-1/4 + a3) + a4 B^-1/4 + a5 )^2
_a1, _a2, _a3, _a4, _a5 = 5.949e-8, -2.389e-7, 2.459e-7, 4.120e-4, -4.225e-4

def R_of_B(B_cd_m2: float) -> float:
    B = float(B_cd_m2)
    if B <= 0:
        B = 1e-12
    term = math.sqrt(_a1 * (B ** (-0.5)) + _a2 * (B ** (-0.25)) + _a3)
    return (term + _a4 * (B ** (-0.25)) + _a5) ** 2


# Full-curve C_infty(B): Eq. 39 with coefficients from your PDF:
# C_infty = sqrt(b1 B^-1/2 + b2 B^-1/4 + b3) + b4 B^-1/4 + b5
_b1, _b2, _b3, _b4, _b5 = 9.606e-6, -6.337e-5, 1.092e-4, 5.791e-3, -6.177e-3

def Cinf_of_B(B_cd_m2: float) -> float:
    B = float(B_cd_m2)
    if B <= 0:
        B = 1e-12
    return math.sqrt(_b1 * (B ** (-0.5)) + _b2 * (B ** (-0.25)) + _b3) + _b4 * (B ** (-0.25)) + _b5


def crumey_threshold_contrast(B_cd_m2: float, A_arcmin2: float, *, field_factor: float = 1.0) -> float:
    """
    Threshold Weber contrast for a target of angular area A (arcmin^2)
    on background luminance B (cd/m^2), using Crumey Eq. 41:
      C = ((R/A)^q + C_infty^q)^(1/q)

    field_factor >= 1 multiplies the threshold to make predictions conservative.
    """
    B = float(B_cd_m2)
    A = float(A_arcmin2)
    if A <= 0:
        raise ValueError("A_arcmin2 must be > 0")
    q = q_of_B(B)
    R = R_of_B(B)
    Cinf = Cinf_of_B(B)
    C = ((R / A) ** q + (Cinf ** q)) ** (1.0 / q)
    return float(field_factor) * max(0.0, C)


# ------------------- Telescope adaptation (Crumey Eq. 66) -------------------
def crumey_Ba(B_cd_m2: float, spec: TelescopeSpec) -> Tuple[float, Dict[str, float]]:
    """
    Effective background luminance at the eye when viewing through a telescope.

    Crumey Eq. 66:
      B_a = (delta_min/p)^2 * B / F_t
    with delta_min = min(d_exit, p) (Eq. 67) and F_t^{-1} = telescope transmittance.

    We use:
      transmittance T = telescope_throughput(spec)
      so B / F_t = B * T
    """
    B = float(B_cd_m2)
    p = eye_pupil_diameter_mm(spec.observer_age)
    d_exit = float(spec.aperture_mm) / float(spec.magnification)
    delta_min = min(d_exit, p)
    T = telescope_throughput(spec)
    Ba = (delta_min / p) ** 2 * B * T
    return Ba, {
        "eye_pupil_mm": p,
        "exit_pupil_mm": d_exit,
        "delta_min_mm": delta_min,
        "throughput_T": T,
    }


def crumey_telescope_threshold_contrast(B_cd_m2: float, A_arcmin2: float, *, spec: TelescopeSpec) -> Tuple[float, Dict[str, float]]:
    """
    Telescope threshold for an extended target:

    Crumey Eq. 77:
      C_tel = phi * ((R_a/A_a)^q + C_a^q)^(1/q)

    In practice we evaluate the same general threshold model at:
      B_a = effective background at eye
      A_a = A * M^2

    and then multiply by:
      spec.field_factor (conservative)
      spec.phi (optional Crumey field factor; set 1.0 if unsure)
    """
    Ba, d0 = crumey_Ba(B_cd_m2, spec)
    Aa = float(A_arcmin2) * (float(spec.magnification) ** 2)

    C_tel = crumey_threshold_contrast(Ba, Aa, field_factor=spec.field_factor) * float(spec.phi)
    d0.update({
        "Ba_cd_m2": Ba,
        "Aa_arcmin2": Aa,
        "q(Ba)": q_of_B(Ba),
        "R(Ba)": R_of_B(Ba),
        "Cinf(Ba)": Cinf_of_B(Ba),
    })
    return C_tel, d0


# ------------------- High-level helper -------------------
def visibility_margin(
    L_nl: float,
    B_nl: float,
    A_arcmin2: float,
    *,
    spec: TelescopeSpec | None = None,
    use_abs_contrast: bool = True,
    naked_eye_field_factor: float = 2.4,
) -> Dict[str, float]:
    """
    Compute Crumey visibility margin (naked-eye + optional telescope).

    margin = log10(C_use / C_thr)
      >0  : detectable (per model)
      <0  : not detectable

    For hilal you usually want use_abs_contrast=True.

    naked_eye_field_factor: multiplier to raise naked-eye threshold from lab data;
    you can tune it to match your local observation statistics.
    """
    if B_nl <= 0:
        raise ValueError("B_nl must be > 0")

    L = float(L_nl)
    B = float(B_nl)
    B_cd = B * NL_TO_CD_M2

    # Weber contrast
    C_obj = (L - B) / B
    C_use = abs(C_obj) if use_abs_contrast else C_obj

    # thresholds
    C_thr_naked = crumey_threshold_contrast(B_cd, A_arcmin2, field_factor=naked_eye_field_factor)
    margin_naked = math.log10(C_use / C_thr_naked) if (C_use > 0 and C_thr_naked > 0) else float("-inf")

    out = {
        "C_obj": C_obj,
        "C_use": C_use,
        "B_cd_m2": B_cd,
        "A_arcmin2": float(A_arcmin2),
        "C_thr_naked": C_thr_naked,
        "margin_naked": margin_naked,
    }

    if spec is not None:
        C_thr_tel, diag = crumey_telescope_threshold_contrast(B_cd, A_arcmin2, spec=spec)
        margin_tel = math.log10(C_use / C_thr_tel) if (C_use > 0 and C_thr_tel > 0) else float("-inf")
        out.update({
            "C_thr_tel": C_thr_tel,
            "margin_tel": margin_tel,
        })
        out.update(diag)

    return out
