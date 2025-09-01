@tool
extends DirectionalLight3D

const HALF_PI_F: float = PI / 2.0

var atmosphere_tex := load(get_script().resource_path.get_base_dir() + "/sky_lut.tres")
var cloud_sky := load(get_script().resource_path.get_base_dir() + "/clouds_sky.tres")
var date_time_location: Node3D

func _ready():
	call_deferred("delayed_init")
	
# Needs to be called after cloud_sky is fully initialized
func delayed_init():
	cloud_sky.sun = self
	cloud_sky.request_full_sky_init()
	date_time_location = cloud_sky.date_time_location

func _process(delta: float) -> void:
	var utc_date_and_time = date_time_location.convert_to_utc_date_and_time()
	var pos = sun_position(utc_date_and_time.year, utc_date_and_time.day_of_the_year, utc_date_and_time.time_in_hours, date_time_location.latitude, date_time_location.longitude)
	rotation.y = HALF_PI_F - pos.azimuth
	rotation.x = -pos.elevation
	atmosphere_tex.request_update()

# Sun position calculations are based on information found in a stack overflow post:
# https://stackoverflow.com/questions/8708048/position-of-the-sun-given-time-of-day-latitude-and-longitude
#
# More sites with information about this:
# https://www.stjarnhimlen.se/comp/tutorial.html#5 sun position calc
# Has a similar calculation and continues with the moon position too, explanations are quite clear.
func sun_position(year: int, day_of_the_year: int, exact_hour: float, lat: float, long: float):
	# Julian date
	var delta = year - 1949
	var leap = floor(delta / 4) # former leapyears
	var julian_date = 32916.5 + delta * 365.0 + leap + day_of_the_year + exact_hour / 24.0

	# The input to the Astronomer's almanach is the difference between
	# the Julian date and JD 2451545.0 (noon, 1 January 2000)
	var time = julian_date - 51545.0

	# Ecliptic coordinates
	# Mean longitude
	var mnlong_deg = fmod(280.460 + 0.9856474 * time, 360.0)

	# Mean anomaly
	var mnanom_rad = deg_to_rad(fmod(357.528 + 0.9856003 * time, 360.0))

	# Ecliptic longitude and obliquity of ecliptic
	var eclong_rad = deg_to_rad(fmod(mnlong_deg + 1.915 * sin(mnanom_rad) + 0.020 * sin(2.0 * mnanom_rad), 360.0))
	var oblqec_rad = deg_to_rad(23.439 - 0.0000004 * time)

	# Celestial coordinates
	# Right ascension and declination
	var num = cos(oblqec_rad) * sin(eclong_rad)
	var den = cos(eclong_rad)
	var ra_rad = atan(num / den)
	if den < 0.0:
		ra_rad = ra_rad + PI
	elif num < 0.0:
		ra_rad = ra_rad + 2.0 * PI
	var dec_rad = asin(sin(oblqec_rad) * sin(eclong_rad))

	# Local coordinates
	# Greenwich mean sidereal time
	var gmst = fmod(6.697375 + 0.0657098242 * time + exact_hour, 24.0)

	# Local mean sidereal time
	var lmst = fmod(gmst + long / 15.0, 24.0)
	var lmst_rad = deg_to_rad(15.0 * lmst)

	# Hour angle (rad)
	var ha_rad = fmod(lmst_rad - ra_rad, (2.0 * PI))

	# Latitude [rad]
	var lat_rad = deg_to_rad(lat)

	# Elevation
	var el_rad = asin(sin(dec_rad) * sin(lat_rad) + cos(dec_rad) * cos(lat_rad) * cos(ha_rad))
	
	# New code: Solar zenith angle
	var zenithAngle = acos(sin(lat_rad) * sin(dec_rad) + cos(lat_rad) * cos(dec_rad) * cos(ha_rad))
	# Solar azimuth
	var az_rad = acos(((sin(lat_rad) * cos(zenithAngle)) - sin(dec_rad)) / (cos(lat_rad) * sin(zenithAngle)))

	# New code
	if (ha_rad > 0):
		az_rad = az_rad + deg_to_rad(180)
	else:
		az_rad = deg_to_rad(540) - az_rad

	# Azimuth
	az_rad = fmod(az_rad, deg_to_rad(360))

	return { "azimuth": az_rad, "elevation": el_rad }
