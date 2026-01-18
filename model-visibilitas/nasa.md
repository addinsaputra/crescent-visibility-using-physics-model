saya lampirkan modul Python **meteo_power_api.py** dan contoh skrip **hilal_with_power_demo.py** yang mengotomatisasi pengambilan data kelembapan relatif dan suhu udara dari **NASA POWER** lalu menginterpolasinya ke waktu senja/sunset yang diinginkan. Modul ini memanfaatkan data *Hourly* dari POWER dengan **time-standard=UTC** (default API adalah LST) agar waktu yang diterima sejalan dengan perhitungan astronomi. Dalam FAQ NASA dijelaskan bahwa cap waktu jamannya mewakili **awal jam**; ini penting untuk interpolasi.

### Struktur data & algoritma

* `ObservingLocation`: kelas data tak berubah yang menyimpan nama lokasi, lintang, bujur, ketinggian dan zona waktu.

* `_fetch_power_hourly(day_utc, location)`: fungsi privat untuk mengirim permintaan ke endpoint `https://power.larc.nasa.gov/api/temporal/hourly/point` dengan parameter `RH2M,T2M`, posisi geografi, rentang tanggal (start=end=YYYYMMDD) dan `time-standard=UTC`. Permintaan dikembalikan dalam format CSV; fungsi ini menangani header metadata sebelum data sebenarnya.

* `_parse_power_csv(...)` dan `_row_to_datetime(...)`: fungsi utilitas untuk mengekstrak baris data dari CSV dan mengubah kolom `YEAR,MO,DY,HR` menjadi objek `datetime` bertimezone UTC.

* `_extract_numeric(...)`: mengonversi nilai string ke float serta mendeteksi kode data hilang `-999`.

* `interpolate_linear(...)`: interpolasi linier sederhana antara dua titik `(x0,y0)` dan `(x1,y1)`.

* `get_rh_t_at_time(location, target_utc)`: fungsi publik yang:

  1. Menentukan dua sampel jam yang mengapit `target_utc` (pembulatan ke bawah jam dan jam berikutnya).
  2. Mengambil data POWER per hari untuk tanggal‐tanggal yang diperlukan (satu permintaan per tanggal).
  3. Membangun kamus `datetime` → baris data.
  4. Mengambil RH2M dan T2M pada jam t₀ dan t₁; jika ada yang hilang, mengeluarkan `PowerAPIError`.
  5. Melakukan interpolasi linier terhadap detik sejak t₀ untuk mendapatkan kelembapan dan suhu di `target_utc`, serta mengapit nilai RH pada 0–100 %.

* Fungsi pengambilan `_fetch_power_hourly` dibungkus dekorator cache sederhana sehingga permintaan untuk kombinasi lokasi dan tanggal yang sama tidak diulang.

Modul ini disertai contoh skrip **hilal_with_power_demo.py**. Skrip ini mendefinisikan lokasi (Semarang, Indonesia), menentukan sunset UTC (misal 10:32 UTC setara 17:32 WIB), memanggil `get_rh_t_at_time` untuk mendapatkan RH dan suhu terinterpolasi, lalu meneruskan nilai tersebut ke `compute_hilal_visibility` di file Anda. Skrip ini menunjukkan cara mengganti input manual kelembapan dan suhu dengan data otomatis.

tolong integrasikan modul ini ke program hilal dalam folder model-visibilitas dengan memanggil `get_rh_t_at_time` sebelum menghitung visibilitas. Pastikan untuk menghitung waktu sunset dalam UTC agar sesuai dengan data POWER. setelah itu uji hingga berhasil.
