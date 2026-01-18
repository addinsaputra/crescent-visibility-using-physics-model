import numpy as np 
from skyfield import api
from skyfield import almanac
from datetime import datetime
from datetime import timedelta
import re
try:
	import pytz
	def _make_tz(tz_in):
		if tz_in is None:
			return pytz.timezone('UTC')
		if isinstance(tz_in, int):
			offset = int(tz_in)
			if offset >= 0:
				tz_name = f"Etc/GMT-{offset}"
			else:
				tz_name = f"Etc/GMT+{abs(offset)}"
			return pytz.timezone(tz_name)
		tzs = str(tz_in).strip()
		if re.fullmatch(r"[+-]?\\d+", tzs):
			offset = int(tzs)
			if offset >= 0:
				tz_name = f"Etc/GMT-{offset}"
			else:
				tz_name = f"Etc/GMT+{abs(offset)}"
			return pytz.timezone(tz_name)
		return pytz.timezone(tzs)

	def _localize(tz_obj, dt):
		return tz_obj.localize(dt)
except Exception:
	from zoneinfo import ZoneInfo

	def _make_tz(tz_in):
		if tz_in is None:
			return ZoneInfo('UTC')
		if isinstance(tz_in, int):
			offset = int(tz_in)
			if offset >= 0:
				tz_name = f"Etc/GMT-{offset}"
			else:
				tz_name = f"Etc/GMT+{abs(offset)}"
			return ZoneInfo(tz_name)
		tzs = str(tz_in).strip()
		if re.fullmatch(r"[+-]?\\d+", tzs):
			offset = int(tzs)
			if offset >= 0:
				tz_name = f"Etc/GMT-{offset}"
			else:
				tz_name = f"Etc/GMT+{abs(offset)}"
			return ZoneInfo(tz_name)
		return ZoneInfo(tzs)

	def _localize(tz_obj, dt):
		return dt.replace(tzinfo=tz_obj)
from skyfield.units import Angle
import os
from open_meteo_api_era5 import fetch_hourly_weather

__all__ = ["list_hijri_months", "hijri_month", "set_location", "convert_utc_to_localtime", "convert_localtime_to_utc", "sunrise_sunset_utc",
			"sunrise_sunset_local", "sunrise_sunset_geometric_local", "sunrise_sunset_apparent_local",
			"sun_position_time_utc", "sun_position_time_local", "moon_position_time_utc", "moon_position_time_local",
			"moon_elongation_time_utc", "moon_elongation_time_local", "moon_phase_angle_time_utc", "moon_phase_angle_time_local",
			"moon_illumination_width_utc", "moon_illumination_width_local",
			"find_new_moon_dates", "ref_hijri_ijtima", "newmoon_hijri_month_utc", "newmoon_hijri_month_local_time", "refraction_horizon_degree",
			"moonrise_moonset_utc", "moonrise_moonset_local", "print_angle", "print_timedelta", "print_timedelta_tz", "fajr_time_utc",
			"fajr_time_local", "calc_timedelta_seconds"]

global ts, ephem
ts = api.load.timescale()
_base_dir = os.path.dirname(__file__)
_ephem_path = os.path.join(_base_dir, 'de421.bsp')
if not os.path.exists(_ephem_path):
	_ephem_path = os.path.join(os.path.dirname(_base_dir), 'de421.bsp')
if not os.path.exists(_ephem_path):
	_ephem_path = os.path.join(os.getcwd(), 'de421.bsp')
ephem = api.load_file(_ephem_path)

def list_hijri_months(print_list=False):
	hijri_months = ['Muharram', 'Shafar', 'Rabiul Awwal', 'Rabiuts Tsani', 'Jumadil Ula', 'Jumadil Akhir', 'Rajab', 'Syaban', 'Ramadhan', 'Syawal', 'Dzulqadah', 'Dzulhijjah']
	
	if print_list == True:
		for ii in range(0,12):
			print ('%d %s' % (ii+1,hijri_months[ii]))

	return hijri_months

def hijri_month(idx):
	hijri_months = list_hijri_months()
	return hijri_months[int(idx)-1]

def set_location(latitude, longitude, elevation):
	if latitude > 0:
		lat_str = '%lf N' % latitude
	else:
		lat_str = '%lf S' % (-1.0*latitude)

	if longitude > 0:
		long_str = '%lf E' % longitude
	else:
		long_str = '%lf W' % (-1.0*longitude)

	location = api.Topos(lat_str, long_str, elevation_m=elevation)
	return location


def get_temperature_from_power_api(skyfield_location, target_utc, timezone_str='UTC', name='Location'):
	"""Get temperature from Open-Meteo ERA5 API at a specific time.

	Parameters
	----------
	skyfield_location : skyfield.toposlib.Topos
		Location object from skyfield set_location function
	target_utc : datetime
		Timezone-aware UTC datetime for which to get temperature
	timezone_str : str, optional
		IANA timezone identifier (default: 'UTC')
	name : str, optional
		Name of the location (default: 'Location')

	Returns
	-------
	float
		Temperature in degrees Celsius

	Raises
	------
	Exception
		If the API request fails or data is missing
	"""
	# Extract coordinates from skyfield location
	latitude = skyfield_location.latitude.degrees
	longitude = skyfield_location.longitude.degrees

	# Ensure target_utc is timezone-aware
	from datetime import timezone as tz
	if target_utc.tzinfo is None:
		target_utc = target_utc.replace(tzinfo=tz.utc)

	# Fetch hourly weather data from ERA5 API
	date_str = target_utc.strftime("%Y-%m-%d")
	df_weather = fetch_hourly_weather(
		latitude=latitude,
		longitude=longitude,
		start_date=date_str,
		end_date=date_str,
		hourly_variables=["temperature_2m"],
		timezone="UTC"
	)

	# Find the data closest to target time
	df_weather['time_diff'] = abs(df_weather['date'] - target_utc)
	idx_closest = df_weather['time_diff'].idxmin()
	temp_c = float(df_weather.loc[idx_closest, 'temperature_2m'])

	return temp_c


def convert_utc_to_localtime(time_zone_str, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	# determine timezone object (supports numeric offsets like "+7" or IANA names)
	time_zone = _make_tz(time_zone_str)

	if utc_datetime is None:
		t = ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)
		utc_datetime = t.utc_datetime()
		local_datetime = utc_datetime.astimezone(time_zone)
	else:
		local_datetime = utc_datetime.astimezone(time_zone)

	return local_datetime


def convert_localtime_to_utc(time_zone_str, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	# determine timezone object (supports numeric offsets like "+7" or IANA names)
	time_zone = _make_tz(time_zone_str)

	if local_datetime is None:
		d = datetime(year, month, day, hour, minute, int(second))
		tz = _localize(time_zone, d)
		t = ts.from_datetime(tz)
		utc_datetime = t.utc_datetime()
	else:
		if local_datetime.tzinfo is None:
			local_datetime = _localize(time_zone, local_datetime)
		t = ts.from_datetime(local_datetime)
		utc_datetime = t.utc_datetime()

	return utc_datetime


def refraction_horizon_degree(temperature_C, pressure_mbar):
	from skyfield.earthlib import refraction

	r = refraction(0.0, temperature_C=temperature_C, pressure_mbar=pressure_mbar)  # in degree
	return r


def sunrise_sunset_geometric_local(location, time_zone_str, year=None, month=None, day=None, radius_degrees=0.2665):
	"""Hitung sunrise/sunset tanpa koreksi refraksi atmosfer (hanya koreksi geometri/elevasi).

	Fungsi ini menghitung sunset berdasarkan koreksi geometri akibat elevasi lokasi,
	TANPA koreksi refraksi atmosfer. Gunakan fungsi ini untuk mendapatkan estimasi
	waktu sunset terlebih dahulu, kemudian ambil data atmosfer, lalu gunakan
	sunrise_sunset_apparent_local() untuk koreksi refraksi.

	Parameters
	----------
	location : skyfield.toposlib.Topos
		Lokasi pengamatan dari set_location()
	time_zone_str : str
		Zona waktu (contoh: "Asia/Jakarta" atau "+7")
	year, month, day : int
		Tanggal untuk perhitungan
	radius_degrees : float
		Radius matahari (default: 0.2665 derajat)

	Returns
	-------
	sunrise_local : datetime
		Waktu sunrise dalam waktu lokal
	sunset_local : datetime
		Waktu sunset dalam waktu lokal (tanpa koreksi refraksi)
	"""
	# get corection due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	# HORIZON GEOMETRIS - TANPA KOREKSI REFRAKSI
	horizon_degrees = h.degrees  # Tidak ada refraksi_horizon_degree

	t0 = ts.utc(year, month, day, 0)
	t1 = ts.utc(t0.utc_datetime() + timedelta(days=2))

	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location,
													horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
	sunrise = None
	sunset = None
	for time, is_sunrise in zip(t, y):
		if is_sunrise:
			sunrise0 = time.utc_datetime()
			if (sunrise0.hour*60.0)+sunrise0.minute+(sunrise0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunrise0.day == day_plus1.day:
					sunrise = sunrise0
			else:
				if sunrise0.day == day:
					sunrise = sunrise0

		else:
			sunset0 = time.utc_datetime()
			if (sunset0.hour*60.0)+sunset0.minute+(sunset0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunset0.day == day_plus1.day:
					sunset = sunset0
			else:
				if sunset0.day == day:
					sunset = sunset0

	# Konversi ke local time
	if sunrise is not None:
		sunrise_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunrise)
	else:
		sunrise_local = None
	sunset_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunset)

	return sunrise_local, sunset_local


def sunrise_sunset_apparent_local(location, time_zone_str, year=None, month=None, day=None,
								   temperature_C=10.0, pressure_mbar=1030.0, radius_degrees=0.2665):
	"""Hitung sunrise/sunset dengan koreksi refraksi atmosfer.

	Fungsi ini menghitung sunset dengan koreksi refraksi atmosfer menggunakan
	algoritma yang sama dengan sunmoon_01.py (almanac.find_discrete).

	ALUR PENGGUNAAN YANG DISARANKAN:
	1. Panggil sunrise_sunset_geometric_local() untuk mendapatkan estimasi sunset
	2. Gunakan waktu tersebut untuk mengambil data atmosfer (T, RH) dari API
	3. Panggil fungsi ini dengan temperature dari API untuk mendapatkan sunset yang akurat

	Parameters
	----------
	location : skyfield.toposlib.Topos
		Lokasi pengamatan dari set_location()
	time_zone_str : str
		Zona waktu (contoh: "Asia/Jakarta" atau "+7")
	year, month, day : int
		Tanggal untuk perhitungan
	temperature_C : float
		Suhu dalam derajat Celsius untuk koreksi refraksi
	pressure_mbar : float
		Tekanan udara dalam mbar untuk koreksi refraksi
	radius_degrees : float
		Radius matahari (default: 0.2665 derajat)

	Returns
	-------
	sunrise_local : datetime
		Waktu sunrise dalam waktu lokal
	sunset_local : datetime
		Waktu sunset dalam waktu lokal dengan koreksi refraksi
	"""
	# get correction due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	# HORIZON DENGAN KOREKSI REFRAKSI (menggunakan algoritma sunmoon_01.py)
	horizon_degrees = h.degrees - refraction_horizon_degree(temperature_C, pressure_mbar)

	t0 = ts.utc(year, month, day, 0)
	t1 = ts.utc(t0.utc_datetime() + timedelta(days=2))

	# Gunakan almanac.find_discrete untuk menemukan sunrise/sunset
	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location,
												horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))

	sunrise = None
	sunset = None
	for time, is_sunrise in zip(t, y):
		if is_sunrise:
			sunrise0 = time.utc_datetime()
			if (sunrise0.hour*60.0)+sunrise0.minute+(sunrise0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunrise0.day == day_plus1.day:
					sunrise = sunrise0
			else:
				if sunrise0.day == day:
					sunrise = sunrise0
		else:
			sunset0 = time.utc_datetime()
			if (sunset0.hour*60.0)+sunset0.minute+(sunset0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunset0.day == day_plus1.day:
					sunset = sunset0
			else:
				if sunset0.day == day:
					sunset = sunset0

	# Konversi ke local time
	if sunrise is not None:
		sunrise_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunrise)
	else:
		sunrise_local = None
	
	if sunset is not None:
		sunset_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunset)
	else:
		sunset_local = None

	return sunrise_local, sunset_local 


def sunrise_sunset_utc(location, year=None, month=None, day=None, temperature_C=10.0, pressure_mbar=1030.0, radius_degrees=0.2665, auto_temperature=False, timezone_str='UTC', location_name='Location'):

	# Get temperature from POWER API if auto_temperature is enabled
	if auto_temperature:
		from datetime import timezone as tz
		# Create timezone-aware datetime for noon UTC on the given date
		target_dt = datetime(year, month, day, 12, 0, 0, tzinfo=tz.utc)
		try:
			temperature_C = get_temperature_from_power_api(location, target_dt, timezone_str, location_name)
		except Exception as e:
			print(f"Warning: Failed to get temperature from POWER API ({e}). Using default value {temperature_C}°C")

	# get corection due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	horizon_degrees = h.degrees - refraction_horizon_degree(temperature_C, pressure_mbar)

	t0 = ts.utc(year, month, day, 0)
	t1 = ts.utc(t0.utc_datetime() + timedelta(days=2))

	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location, 
													horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
	#print (t.utc_datetime())
	sunrise = None
	sunset = None
	for time, is_sunrise in zip(t, y):
		if is_sunrise:
			sunrise0 = time.utc_datetime()
			if (sunrise0.hour*60.0)+sunrise0.minute+(sunrise0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunrise0.day == day_plus1.day:
					sunrise = sunrise0
			else:
				if sunrise0.day == day:
					sunrise = sunrise0

		else:
			sunset0 = time.utc_datetime()
			if (sunset0.hour*60.0)+sunset0.minute+(sunset0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if sunset0.day == day_plus1.day:
					sunset = sunset0
			else:
				if sunset0.day == day:
					sunset = sunset0

	return sunrise, sunset


def sunrise_sunset_local(location, time_zone_str, year=None, month=None, day=None, temperature_C=10.0, pressure_mbar=1030.0, radius_degrees=0.2665, auto_temperature=False, location_name='Location'):
	sunrise_utc, sunset_utc = sunrise_sunset_utc(location, year=year, month=month, day=day,
												temperature_C=temperature_C, pressure_mbar=pressure_mbar, radius_degrees=radius_degrees,
												auto_temperature=auto_temperature, timezone_str=time_zone_str, location_name=location_name)
	sunrise_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunrise_utc)
	sunset_local = convert_utc_to_localtime(time_zone_str, utc_datetime=sunset_utc)

	return sunrise_local, sunset_local


def fajr_time_utc(location, year=None, month=None, day=None, temperature_C=10.0, pressure_mbar=1030.0, fajr_sun_altitude=-18.0):

	# get corection due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	horizon_degrees = h.degrees - refraction_horizon_degree(temperature_C, pressure_mbar)

	radius_degrees = horizon_degrees -1*fajr_sun_altitude

	t0 = ts.utc(year, month, day, 0)
	t1 = ts.utc(t0.utc_datetime() + timedelta(days=1))

	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location, horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
	fajr_time = None
	for time, is_sunrise in zip(t, y):
		if is_sunrise:
			fajr_time = time.utc_datetime()

	if fajr_time is None:
		t0 = ts.utc(year, month, day, 0)
		t1 = ts.utc(t0.utc_datetime() + timedelta(days=2))

		t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location, horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
		for time, is_sunrise in zip(t, y):
			if is_sunrise:
				fajr_time = time.utc_datetime()

	return fajr_time


def fajr_time_local(location, time_zone_str, year=None, month=None, day=None, temperature_C=10.0, pressure_mbar=1030.0, fajr_sun_altitude=-18.0):

	# get corection due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	horizon_degrees = h.degrees - refraction_horizon_degree(temperature_C, pressure_mbar)

	#t0 = ts.utc(year, month, day, 0)
	#t1 = ts.utc(t0.utc_datetime() + timedelta(days=1))

	radius_degrees = horizon_degrees -1*fajr_sun_altitude

	tinit = ts.utc(year, month, day, 0)
	t0 = ts.utc(tinit.utc_datetime() - timedelta(days=1))
	t1 = ts.utc(tinit.utc_datetime() + timedelta(days=1))

	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Sun'], location, horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
	fajr_time = None
	for time, is_sunrise in zip(t, y):
		if is_sunrise:
			fajr_utc = time.utc_datetime()
			fajr_local = convert_utc_to_localtime(time_zone_str, utc_datetime=fajr_utc)
			if fajr_local.day == day:
				fajr_time = fajr_local

	return fajr_time


def moonrise_moonset_utc(location, year=None, month=None, day=None, temperature_C=23.0, pressure_mbar=1030.0, radius_degrees=0.2575, auto_temperature=False, timezone_str='UTC', location_name='Location'):

	# Get temperature from POWER API if auto_temperature is enabled
	if auto_temperature:
		from datetime import timezone as tz
		# Create timezone-aware datetime for noon UTC on the given date
		target_dt = datetime(year, month, day, 12, 0, 0, tzinfo=tz.utc)
		try:
			temperature_C = get_temperature_from_power_api(location, target_dt, timezone_str, location_name)
		except Exception as e:
			print(f"Warning: Failed to get temperature from POWER API ({e}). Using default value {temperature_C}°C")

	# get corection due to elevation
	elev = location.elevation
	altitude_m = elev.m
	earth_radius_m = 6378136.6
	side_over_hypotenuse = earth_radius_m/(earth_radius_m + altitude_m)
	h = Angle(radians = -np.arccos(side_over_hypotenuse))

	horizon_degrees = h.degrees - refraction_horizon_degree(temperature_C, pressure_mbar)

	tinit = ts.utc(year, month, day, 0)
	t0 = ts.utc(tinit.utc_datetime() - timedelta(days=1))
	t1 = ts.utc(tinit.utc_datetime() + timedelta(days=3))

	t, y = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephem, ephem['Moon'], location, 
													horizon_degrees=horizon_degrees, radius_degrees=radius_degrees))
	#print (t.utc_datetime())
	moonrise = None
	moonset = None
	for time, is_moonrise in zip(t, y):
		if is_moonrise:
			moonrise0 = time.utc_datetime()
			if (moonrise0.hour*60.0)+moonrise0.minute+(moonrise0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if moonrise0.day == day_plus1.day:
					moonrise = moonrise0
			else:
				if moonrise0.day == day:
					moonrise = moonrise0

		else:
			moonset0 = time.utc_datetime()
			if (moonset0.hour*60.0)+moonset0.minute+(moonset0.second/60.0) + location.longitude.degrees*4 < 0.0:
				day_plus1 = datetime(year,month,day,0,0,0) + timedelta(days=1)
				if moonset0.day == day_plus1.day:
					moonset = moonset0
			else:
				if moonset0.day == day:
					moonset = moonset0

	return moonrise, moonset


def moonrise_moonset_local(location, time_zone_str, year=None, month=None, day=None, temperature_C=22.0, pressure_mbar=1030.0, radius_degrees=0.2575, auto_temperature=False, location_name='Location'):
	moonrise_utc, moonset_utc = moonrise_moonset_utc(location, year=year, month=month, day=day,
									temperature_C=temperature_C, pressure_mbar=pressure_mbar, radius_degrees=radius_degrees,
									auto_temperature=auto_temperature, timezone_str=time_zone_str, location_name=location_name)
	if moonrise_utc is not None:
		moonrise_local = convert_utc_to_localtime(time_zone_str, utc_datetime=moonrise_utc)
	else:
		moonrise_local = None
	moonset_local = convert_utc_to_localtime(time_zone_str, utc_datetime=moonset_utc)

	return moonrise_local, moonset_local


def sun_position_time_utc(location, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None, temperature_C=10.0, pressure_mbar=1030.0):
	sun = ephem["Sun"]
	earth = ephem["Earth"]

	if utc_datetime is None:
		sun_pos = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(sun).apparent()
	else:
		sun_pos = (earth + location).at(ts.from_datetime(utc_datetime)).observe(sun).apparent()

	altitude, azimuth, distance = sun_pos.altaz(temperature_C=temperature_C, pressure_mbar=pressure_mbar)

	return altitude.degrees, azimuth.degrees, distance.km


def sun_position_time_local(location, time_zone_str, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None, temperature_C=10.0, pressure_mbar=1030.0):

	utc_datetime = convert_localtime_to_utc(time_zone_str, local_datetime=local_datetime, 
						year=year, month=month, day=day, hour=hour, minute=minute, second=second)

	sun_alt, sun_az, distance = sun_position_time_utc(location, utc_datetime=utc_datetime, temperature_C=temperature_C, pressure_mbar=pressure_mbar)

	return sun_alt, sun_az, distance


def moon_position_time_utc(location, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None, temperature_C=10.0, pressure_mbar=1030.0):
	moon = ephem["Moon"]
	earth = ephem["Earth"]

	if utc_datetime is None:
		moon_pos = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(moon).apparent()
	else:
		moon_pos = (earth + location).at(ts.from_datetime(utc_datetime)).observe(moon).apparent()

	altitude, azimuth, distance = moon_pos.altaz(temperature_C=temperature_C, pressure_mbar=pressure_mbar)

	return altitude.degrees, azimuth.degrees, distance.km


def moon_position_time_local(location, time_zone_str, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None, temperature_C=10.0, pressure_mbar=1030.0):

	utc_datetime = convert_localtime_to_utc(time_zone_str, local_datetime=local_datetime, 
						year=year, month=month, day=day, hour=hour, minute=minute, second=second)

	moon_alt, moon_az, distance = moon_position_time_utc(location, utc_datetime=utc_datetime, temperature_C=temperature_C, pressure_mbar=pressure_mbar)

	return moon_alt, moon_az, distance


def moon_elongation_time_utc(location=None, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	earth = ephem["Earth"]
	moon = ephem["Moon"]
	sun = ephem["Sun"]

	if utc_datetime is None:
		if location is None:
			m = earth.at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(moon).apparent()
			s = earth.at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(sun).apparent()
		else:
			m = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(moon).apparent()
			s = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(sun).apparent()
	else:
		if location is None:
			m = earth.at(ts.from_datetime(utc_datetime)).observe(moon).apparent()
			s = earth.at(ts.from_datetime(utc_datetime)).observe(sun).apparent()
		else:
			m = (earth + location).at(ts.from_datetime(utc_datetime)).observe(moon).apparent()
			s = (earth + location).at(ts.from_datetime(utc_datetime)).observe(sun).apparent()

	return s.separation_from(m).degrees


def moon_elongation_time_local(time_zone_str, location=None, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	
	utc_datetime = convert_localtime_to_utc(time_zone_str, local_datetime=local_datetime, 
						year=year, month=month, day=day, hour=hour, minute=minute, second=second)

	moon_elong = moon_elongation_time_utc(location=location, utc_datetime=utc_datetime)
	return moon_elong


def moon_phase_angle_time_utc(location=None, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	earth = ephem["Earth"]
	moon = ephem["Moon"]
	sun = ephem["Sun"]

	if utc_datetime is None:
		t = ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)
	else:
		t = ts.from_datetime(utc_datetime)

	if location is None:
		observer = earth
	else:
		observer = earth + location

	moon_pos = moon.at(t).position.km
	sun_pos = sun.at(t).position.km
	observer_pos = observer.at(t).position.km

	vec_moon_to_sun = sun_pos - moon_pos
	vec_moon_to_obs = observer_pos - moon_pos

	dot = np.einsum('i...,i...', vec_moon_to_sun, vec_moon_to_obs)
	norm_sun = np.einsum('i...,i...', vec_moon_to_sun, vec_moon_to_sun)
	norm_obs = np.einsum('i...,i...', vec_moon_to_obs, vec_moon_to_obs)
	cos_angle = dot / np.sqrt(norm_sun * norm_obs)
	cos_angle = np.clip(cos_angle, -1.0, 1.0)
	return np.degrees(np.arccos(cos_angle))


def moon_phase_angle_time_local(time_zone_str, location=None, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):
	
	utc_datetime = convert_localtime_to_utc(time_zone_str, local_datetime=local_datetime, 
						year=year, month=month, day=day, hour=hour, minute=minute, second=second)

	phase_angle = moon_phase_angle_time_utc(location=location, utc_datetime=utc_datetime)
	return phase_angle


def moon_illumination_width_utc(location=None, utc_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):

	from pymeeus.Epoch import Epoch
	from pymeeus.Moon import Moon

	earth = ephem["Earth"]
	moon = ephem["Moon"]
	sun = ephem["Sun"]

	if utc_datetime is None:
		if location is None:
			m = earth.at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(moon).apparent()
			s = earth.at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(sun).apparent()
		else:
			m = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(moon).apparent()
			s = (earth + location).at(ts.utc(year, month=month, day=day, hour=hour, minute=minute, second=second)).observe(sun).apparent()
	else:
		if location is None:
			m = earth.at(ts.from_datetime(utc_datetime)).observe(moon).apparent()
			s = earth.at(ts.from_datetime(utc_datetime)).observe(sun).apparent()
		else:
			m = (earth + location).at(ts.from_datetime(utc_datetime)).observe(moon).apparent()
			s = (earth + location).at(ts.from_datetime(utc_datetime)).observe(sun).apparent()

		year, month, day = utc_datetime.year, utc_datetime.month, utc_datetime.day 

	elongation = s.separation_from(m).degrees
	illumination = m.fraction_illuminated(sun)*100   # in percent

	# get horizontal parallax
	epoch = Epoch(year, month, day)
	Lambda, Beta, Delta, parallax = Moon.apparent_ecliptical_pos(epoch)  # parallax in degree
	SD = 0.27254*parallax*60.0    # in arcmin
	width = SD*(1.0 - np.cos(elongation*np.pi/180.0))/60.0   # in degree

	return illumination, width, parallax, SD/60.0


def moon_illumination_width_local(time_zone_str, location=None, local_datetime=None, year=None, month=None, day=None, hour=None, minute=None, second=None):

	utc_datetime = convert_localtime_to_utc(time_zone_str, local_datetime=local_datetime, 
						year=year, month=month, day=day, hour=hour, minute=minute, second=second)

	illumination, width, parallax, SD = moon_illumination_width_utc(location=location, utc_datetime=utc_datetime)

	return illumination, width, parallax, SD


def find_new_moon_dates(start_year, start_month, start_day, end_year, end_month, end_day):
	t0 = ts.utc(start_year, start_month, start_day)
	t1 = ts.utc(end_year, end_month, end_day)

	t, y = almanac.find_discrete(t0, t1, almanac.moon_phases(ephem))
	idx_nm = np.where(y == 0)
	n_nm = len(idx_nm[0])

	new_moon_datetime = []
	if n_nm>0:
		datetime = t.utc_datetime()
		for ii in range(n_nm):
			new_moon_datetime.append(datetime[idx_nm[0][ii]])
	else:
		print ("No New Moon phase is found within the given time range...")

	return new_moon_datetime

def ref_hijri_ijtima():
	# 2022-07-28 --> Muharram 1444
	hijri_m, hijri_y = 1, 1444
	utc_datetime = find_new_moon_dates(2022, 7, 28, 2022, 7, 30)
	return hijri_m, hijri_y, utc_datetime[0]

def newmoon_hijri_month_utc(hijri_year, hijri_month):
	""" Function to find the date of new moon associated with given Hijri month and year
	:param hijri_year:

	:param hijri_month:
		Month number start from 1.
	"""
	ref_hijri_m, ref_hijri_y, ref_utc_datetime = ref_hijri_ijtima()

	if hijri_year > ref_hijri_y:
		stat_forward = 1
	elif hijri_year < ref_hijri_y:
		stat_forward = -1
	else:
		if hijri_month > ref_hijri_m:
			stat_forward = 1
		elif hijri_month < ref_hijri_m:
			stat_forward = -1
		else:
			stat_forward = 0

	if stat_forward == 0:
		return ref_utc_datetime
	elif stat_forward == 1:
		day_minus1 = datetime(ref_utc_datetime.year,ref_utc_datetime.month,ref_utc_datetime.day,0,0,0) - timedelta(days=1)
		new_moon_datetimes = find_new_moon_dates(day_minus1.year, day_minus1.month, day_minus1.day, ref_utc_datetime.year+2*(hijri_year-ref_hijri_y+2), 12, 30)

		idx = np.arange(0,len(new_moon_datetimes))
		idx1 = np.where((idx%12==hijri_month-1) & ((idx/12).astype(int)==hijri_year-ref_hijri_y))
		utc_datetime = new_moon_datetimes[idx1[0][0]]
		return utc_datetime

	else:
		day_plus1 = datetime(ref_utc_datetime.year,ref_utc_datetime.month,ref_utc_datetime.day,0,0,0) + timedelta(days=1)
		new_moon_datetimes = find_new_moon_dates(ref_utc_datetime.year-2*(ref_hijri_y-hijri_year+2), 1, 1, day_plus1.year, day_plus1.month, day_plus1.day)

		idx = np.arange(1,len(new_moon_datetimes)+1) - len(new_moon_datetimes)
		idx1 = np.where((idx%12==hijri_month-1) & ((idx/12).astype(int)-1==hijri_year-ref_hijri_y))
		utc_datetime = new_moon_datetimes[idx1[0][0]]
		return utc_datetime

def newmoon_hijri_month_local_time(hijri_year, hijri_month, time_zone_str):
	utc_datetime = newmoon_hijri_month_utc(hijri_year, hijri_month)
	local_time = convert_utc_to_localtime(time_zone_str, utc_datetime=utc_datetime)
	return local_time

def calc_timedelta_seconds(datetime1, datetime2):
	timedelta = datetime2 - datetime1
	timedelta_s = (timedelta.days*86400.0) + timedelta.seconds
	return timedelta_s

def print_angle(angle_degree):
	return Angle(degrees=angle_degree).dstr(format=u'{0}{1}°{2:02}′{3:02}.{4:0{5}}″')

def print_timedelta(s):
    hours = s // 3600 
    s = s - (hours * 3600)
    minutes = s // 60
    seconds = s - (minutes * 60)
    return '{:02}:{:02}:{:02}'.format(int(hours), int(minutes), int(seconds))

def print_timedelta_tz(s): 
	hours = s // 3600 
	s = s - (hours * 3600)
	minutes = s // 60

	return '{:02}:{:02}'.format(int(hours), int(minutes))



