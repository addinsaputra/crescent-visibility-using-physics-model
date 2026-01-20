# =============================================================================
# DATA LOKASI PENGAMATAN HILAL
# =============================================================================
# Format: Nama Lokasi, Latitude (Lintang), Longitude (Bujur), Elevasi (mdpl)
# 
# Catatan:
# - Latitude positif = Utara, negatif = Selatan
# - Longitude positif = Timur, negatif = Barat
# - Elevasi dalam meter di atas permukaan laut (mdpl)
# =============================================================================

import io
import csv

# Data lokasi dalam format CSV string
DATA_LOKASI_CSV = """
Nabire, -3.36583, 135.49383, 20
Biak, -1.18567, 136.08717, 20
Merauke, -8.49383, 140.40067, 20
POB Lhok Nga, 5.46667, 95.24233, 20
Tower Hilal Lhoong, 5.18367, 95.29850, 20
Pantai Anyer, -6.67033, 105.88500, 20
POB Cibeas Pel. Ratu, -7.07400, 106.53133, 20
Tower Hilal Cikelet, -7.59367, 107.62350, 20
POB Syekh Bela-Belu, -7.73983, 110.35017, 20
Pantai Ngliyep Malang, -8.35000, 112.43333, 20
POB Condrodipo, -7.16967, 112.61733, 20
Pero Konda, -9.60683, 118.98533, 20
Tower Hilal Sulamu, -10.14333, 123.62667, 20
Tower Hilal Meras, 1.48033, 124.83367, 20
Tower Hilal Ternate, -0.79983, 127.29483, 20
Pantai Lhoknga-Aceh Besar, 5.46677, 95.24222, 20
Rooftop Observatorium UIN WS, -6.99167, 110.34806, 20
Menara Hilal desa Marana Donggala-PALU, -0.57861, 119.79070, 20
Gedung BMKG NTT-Kupang, -10.15278, 123.60833, 20
Mega Trade Center-Manado, 1.48028, 124.85000, 20
Pantai Loang Baloq-MATARAM, -8.60408, 116.07440, 20
"""


def get_list_lokasi():
    """
    Mengambil daftar lokasi pengamatan sebagai list of dictionaries.
    
    Returns:
    --------
    list[dict]
        List berisi dictionary dengan keys: nama, lat, lon, elevasi
        Contoh:
        [
            {"nama": "Nabire", "lat": -3.36583, "lon": 135.49383, "elevasi": 20},
            ...
        ]
    """
    f = io.StringIO(DATA_LOKASI_CSV.strip())
    reader = csv.reader(f, skipinitialspace=True)
    
    list_lokasi = []
    for row in reader:
        if len(row) >= 4:
            data = {
                "nama": row[0].strip(),
                "lat": float(row[1]),
                "lon": float(row[2]),
                "elevasi": int(row[3])
            }
            list_lokasi.append(data)
    
    return list_lokasi


def get_lokasi_by_index(index: int):
    """
    Mengambil lokasi berdasarkan index (1-indexed untuk user-friendly).
    
    Parameters:
    -----------
    index : int
        Index lokasi (dimulai dari 1)
        
    Returns:
    --------
    dict or None
        Dictionary lokasi atau None jika index tidak valid
    """
    list_lokasi = get_list_lokasi()
    if 1 <= index <= len(list_lokasi):
        return list_lokasi[index - 1]
    return None


def get_lokasi_by_name(nama: str):
    """
    Mencari lokasi berdasarkan nama (partial match, case-insensitive).
    
    Parameters:
    -----------
    nama : str
        Nama lokasi yang dicari
        
    Returns:
    --------
    dict or None
        Dictionary lokasi pertama yang cocok atau None jika tidak ditemukan
    """
    list_lokasi = get_list_lokasi()
    nama_lower = nama.lower()
    
    for loc in list_lokasi:
        if nama_lower in loc["nama"].lower():
            return loc
    return None


def print_daftar_lokasi():
    """
    Menampilkan daftar lokasi dalam format tabel yang rapi.
    """
    list_lokasi = get_list_lokasi()
    
    print("\n" + "=" * 80)
    print("DAFTAR LOKASI PENGAMATAN HILAL")
    print("=" * 80)
    print(f"{'NO':<4} | {'NAMA LOKASI':<42} | {'LAT':<10} | {'LON':<10} | {'ELV':<5}")
    print("-" * 80)
    
    for i, data in enumerate(list_lokasi, 1):
        print(f"{i:<4} | {data['nama']:<42} | {data['lat']:<10.5f} | {data['lon']:<10.5f} | {data['elevasi']}")
    
    print("=" * 80)
    print(f"Total: {len(list_lokasi)} lokasi")
    return list_lokasi


def pilih_lokasi_interaktif():
    """
    Menampilkan menu interaktif untuk memilih lokasi.
    
    Returns:
    --------
    dict or None
        Dictionary lokasi yang dipilih atau None jika dibatalkan
    """
    list_lokasi = print_daftar_lokasi()
    
    print("\nMasukkan nomor lokasi (1-{}) atau 0 untuk input manual:".format(len(list_lokasi)))
    
    try:
        pilihan = int(input("Pilihan Anda: "))
        
        if pilihan == 0:
            # Input manual
            print("\n--- INPUT LOKASI MANUAL ---")
            nama = input("Nama lokasi: ")
            lat = float(input("Latitude (contoh: -6.9167): "))
            lon = float(input("Longitude (contoh: 110.3480): "))
            elv = int(input("Elevasi (meter): "))
            
            return {
                "nama": nama,
                "lat": lat,
                "lon": lon,
                "elevasi": elv
            }
        
        elif 1 <= pilihan <= len(list_lokasi):
            lokasi = list_lokasi[pilihan - 1]
            print(f"\n✓ Lokasi dipilih: {lokasi['nama']}")
            return lokasi
        
        else:
            print("[!] Nomor tidak valid!")
            return None
            
    except ValueError:
        print("[!] Input harus berupa angka!")
        return None


# --- Testing jika dijalankan langsung ---
if __name__ == "__main__":
    # Test: tampilkan daftar
    print_daftar_lokasi()
    
    # Test: ambil by index
    print("\n--- Test get by index ---")
    lok = get_lokasi_by_index(1)
    print(f"Index 1: {lok}")
    
    # Test: cari by name
    print("\n--- Test get by name ---")
    lok = get_lokasi_by_name("UIN")
    print(f"Cari 'UIN': {lok}")