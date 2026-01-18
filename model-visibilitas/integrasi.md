BUATKAN INTEGRASI PROGRAM DAN BUAT FILE BARU DARI PERINTAH DIBAWAH INI:

ALUR PROGRAM MODEL VISIBILITAS HILAL MODEL KASTNER YANG DIMODIFIKASI:
disclaimer : 
1. model schaefer sebenarnya bukan model yang dikhususkan untuk objek hilal, tapi untuk benda langit yang berada di dekat ufuk seperti planet, satelit, komet, dll. nah dalam model tersebut mempertimbangkan input altitude/azimuth matahari dan bulan dan altitude benda langit tujuan untuk mendapatkan nilai sky brightness benda langit tujuan. implementasi model tersebut untuk kasus hilal akan menimbulkan masalah, karena benda langit tujuan adalah bulan (hilal) itu sendiri. jika input altitude bulan disamakan dengan input altitude benda langit tujuan yang adalah hilal itu sendiri maka akan menimbulkan error nilai sky brightness tak terhingga. begitupun juga sama untuk input selisih azimuth bulan dan azimuth benda langit tujuan yang juga digunakan untuk menentukan nilai sky brighness. agar model tersebut dapat digunakan dengan akurat untuk kasus hilal, maka input nilai altitude moon dikunci nilainya di angka -90 derajat, sedangkan selisih azimuth bulan dan azimuth benda langit tujuan dikunci nilainya di angka 180 derajat.
2. sebenarnya saya sudah melakukan integrasi model dalam file  @finall_rumus.py, tapi input data matahari-bulan dan atmosfer di program tersebut masih manual, nah tugas disini adalah membuat program model visibilitas hilal dengan input data yang dimuat secara otomatis. 


RANCANGAN INTEGRASI PROGRAM MODEL VISIBILITAS HILAL KASTNER YANG DIMODIFIKASI

- menentukan awal bulan hijriah apa yang ingin diprediksi visibilitas hilal-nya (contoh: Muharram 1444) (input user)
- menentukan lokasi (lintang, bujur, elevasi, dan time zone) (input user)
- menghitung waktu ijtima atau konjungsi terjadi kapan saat penentuan awal masuk bulan muharram 1444 (contoh: 29 juli 2022 jam 00:55:5 LT dan 28 juli 2022 jam 17:36:2 UTC) (otomatisasi)
- menghitung waktu sunset (jika waktu ijtima > sunset (contoh ijtima tanggal 28 juli jam 19:30, sunset 17:38), maka perhitungan prediksi visibilitas hilal dilakukan ditanggal berikutnya yaitu 29 juli. jika waktu ijtima < sunset (contoh ijtima tanggal 19 februari jam 03:43, sunset jam 17:43), maka perhitungan prediksi visibilitas hilal dilakukan pada hari ijtima tersebut yaitu 19 februari) (otomatisasi)
- selanjutnya menghitung input data yang diperlukan dalam model prediksi visibilitas ini saat sunset. berikut input data yang diperlukan:

INPUT DALAM MODEL SCHAEFER program berada di file @visual_limit_calc2.py YAITU: 
- Sun altitude (otomatisasi),
- moon altitude (tidak perlu diotomatisasi karena nilainya DIKUNCI di angka -90 derajat) (statis), 
- selisih azimuth bulan dan azimuth benda langit tujuan (untuk mengatasi error karena benda langit tujuan adalah bulan itu sendiri, maka nilainya tidak perlu diotomatisasi, nilainya DIKUNCI di angka 180 derajat untuk menghindari error) (statis), 
- selisih azimuth matahari dan azimuth benda langit tujuan (karena benda langit tujuan adalah bulan itu sendiri, maka selisih azimuth matahari dan azimuth benda langit tujuan adalah selisih azimuth matahari dan azimuth bulan) (otomatisasi),
- altitude benda langit tujuan (karena benda langit adalah hilal itu sendiri, maka altitude benda langit tujuan adalah altitude hilal) (otomatisasi),
- latitude (lintang) tempat (input user),
- elevasi tempat (meter diatas permukaan laut) (input user),
- relative humadity (rh) (otomatisasi_1),
- temperature (t) (otomatisasi_1),
- snellen ratio (nilainya dikunci di angka 1) (statis),
- bulan dan tahun (input user),
- phase of the moon (new, crescent, quarter, gibbous, full), disini nilainya dikunci diangka 0 (statis),

DATA INPUT DARI MODEL KASTNER program berada di file @kastner.py YAITU:
- Moon phase angle (yaitu nilai fase bulan yang dihitung dalam derajat) (otomatisasi),
- elongation toposentrik (otomatisasi),
- Moon semidiameter (otomatisasi), 
- Zenith distance (dihitung dengan cara 90 - altitude benda langit tujuan) (otomatisasi),

* otomatisasi -> nilai input akan dihandle oleh program secara otomatis
* otomatisasi_1 -> nilai input akan dihandle oleh program secara otomatis
* input user -> nilai input diisi secara manual oleh user
* statis -> nilai input jangan diubah-ubah, dibiarkan default

- nilai input otomatisasi input ijtima, sunset, matahari dan bulan dihitung menggunakan program dalam file @sunmoon.py
- nilai input otomatisasi_1 input atmosfer RH dan T, dihitung menggunakann program dalam file @meteo_power_api.py

Menghitung model prediksi visibilitas hilal saat sunset:
Model modifikasi visibilitas hilal modus penglihatan naked eye dihitung dengan rumus berikut:
- Delta_m = 2.5 log(R)
- dengan R adalah Rasio Kontras, dihitung dengan cara Luminansi Hilal di dalam atmosfer dalam satuan nL (L) dibagi sky brightness dalam satuan nL (B) -> (L/B)
- parameter sky brightness nL didapatkan dari program file @visual_limit_calc2.py
- parameter Luminansi Hilal nL didapatkan dari program file @kastner.py

Selanjutnya model modifikasi visibilitas hilal dikoreksi untuk modus penglihatan menggunakan alat teleskop. mekanismenya, nilai sky brightness nL (B) dan luminansi hilal nL (L) masing-masing dikoreksi menggunakan rumus pada program @limit_teleskop. selanjutnya nilai yang sudah dikoreksi teleskop tersebut dihitung dengan rumus Delta_m = 2.5 log(R)
- dengan R adalah Rasio Kontras, L koreksi teleskop dibagi B koreksi teleskop (L*/B*)

catatan : pada program file visual_limit_calc2.py menghitung banyak parameter, tapi yang diambil adalah sky_brightness nL dan Koefisien ekstingsi band V atau k_v (array 3), parameter k_v berperan sangat krusial dalam perhitungan visibilitas hilal kastner modifikasi ini karena digunakan untuk menentukan nilai liminansi hilal nL.

setelah selesai uji coba dengan input berikut:
lokasi  : UIN Walisongo
lintang : -7
bujur   : 110.400000
elevasi : 89
bulan   : Muharram
tahun   : 1444
delta day offset : 0


output yang diharapkan dari program ini adalah variabel-variabel dibawah ini:
- nama tempat
- koordinat tempat (lintang, bujur, dan elevasi)
- time zone
- ijitma terjadi pada tanggal (hari, bulan, tahun masehi) dan jam... (local time dan UTC)
- sunset (matahari terbenam) pada jam...
- delta day offsite
- moon altitude
- elongation toposentrik
- relative humadity (rh)
- temperature (t)
- Luminansi Hilal (nL) Naked Eye
- Sky Brightness (nL) Naked Eye
- Luminansi Hilal (nL) Teleskop  
- Sky Brightness (nL) Teleskop
- Koefisien Extingsi (k_V)
- Rasio Kontras R Naked Eye
- Rasio Kontras R Teleskop
- Delta m (Naked Eye)
- Delta m (Teleskop)
- berikan informasi jika delta m kurang dari 0 maka hilal sulit untuk dilihat