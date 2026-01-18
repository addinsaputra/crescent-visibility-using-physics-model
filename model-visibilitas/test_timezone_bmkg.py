"""Test timezone BMKG API"""
from bmkg_api_forecast import _fetch_bmkg_forecast

# Test Kupang (WIT +9)
print('=' * 60)
print('Testing KUPANG (WIT - UTC+9)')
print('=' * 60)

try:
    data = _fetch_bmkg_forecast('53.71.06.1002')  # Kupang
    cuaca = data['data'][0]['cuaca'][0][:3]
    for rec in cuaca:
        utc = rec['utc_datetime']
        local = rec['local_datetime']
        print(f"  UTC: {utc}  |  Local: {local}")
except Exception as e:
    print(f'Error: {e}')

print()

# Test Semarang (WIB +7)
print('=' * 60)
print('Testing SEMARANG (WIB - UTC+7)')
print('=' * 60)

try:
    data = _fetch_bmkg_forecast('33.74.10.1003')  # Semarang
    cuaca = data['data'][0]['cuaca'][0][:3]
    for rec in cuaca:
        utc = rec['utc_datetime']
        local = rec['local_datetime']
        print(f"  UTC: {utc}  |  Local: {local}")
except Exception as e:
    print(f'Error: {e}')

print()
print("KESIMPULAN:")
print("- Jika selisih Local-UTC = 9 jam untuk Kupang -> WIT benar")
print("- Jika selisih Local-UTC = 7 jam untuk Semarang -> WIB benar")
