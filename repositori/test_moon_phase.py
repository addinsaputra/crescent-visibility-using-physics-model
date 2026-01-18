
from datetime import datetime, timezone
import sys
import os

# Add current dir to path
sys.path.append(os.getcwd())

from vislimit_reconstruction import _moon_phase_angle_deg_from_datetime, _load_skyfield

def test():
    print("--- Testing Skyfield Loading ---")
    sf = _load_skyfield()
    if sf:
        print("Skyfield loaded successfully.")
    else:
        print("Skyfield NOT loaded (fallback mode).")

    print("\n--- Testing Phase Calculation ---")
    # Test a known New Moon approx: Jan 6 2000
    dt_new = datetime(2000, 1, 6, 18, 14, tzinfo=timezone.utc)
    phase_new = _moon_phase_angle_deg_from_datetime(dt_new)
    print(f"Phase at New Moon (approx): {phase_new:.4f} (Expected ~0 or ~360)")

    # Test a known Full Moon approx: Jan 21 2000
    dt_full = datetime(2000, 1, 21, 4, 40, tzinfo=timezone.utc)
    phase_full = _moon_phase_angle_deg_from_datetime(dt_full)
    print(f"Phase at Full Moon (approx): {phase_full:.4f} (Expected ~180)")

    print("\n--- Testing Cache ---")
    sf1 = _load_skyfield()
    sf2 = _load_skyfield()
    if sf1 and sf2 and sf1 is sf2:
        print("Cache IS working (objects are identical).")
    elif sf1 is None:
        print("Cache skipped (Skyfield not loaded).")
    else:
        print("Cache FAILED (objects are different).")

if __name__ == "__main__":
    test()
