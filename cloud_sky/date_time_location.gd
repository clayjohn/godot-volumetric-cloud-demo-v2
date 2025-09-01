## Date and time node that can be used to both set a specific date, time and location and to advance time at a configured speed.
@tool extends Node3D

const SECONDS_PER_MINUTE_F: float = 60.0
const SECONDS_PER_MINUTE_I: int = int(SECONDS_PER_MINUTE_F)

const SECONDS_PER_HOUR_F: float = 60.0 * SECONDS_PER_MINUTE_F
const SECONDS_PER_HOUR_I: int = int(SECONDS_PER_HOUR_F)

const MINUTES_PER_HOUR_F: float = 60.0
const MINUTES_PER_HOUR_I: int = int(MINUTES_PER_HOUR_F)

const SECONDS_PER_DAY_F: float = 24.0 * SECONDS_PER_HOUR_F
const SECONDS_PER_DAY_I: int = int(SECONDS_PER_DAY_F)

const SECONDS_TO_HOURS_F: float = 1.0 / 3600.0

const CUMULATIVE_DAYS_IN_EACH_MONTH : Array[int] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
const CUMULATIVE_DAYS_IN_EACH_MONTH_IN_LEAP_YEAR : Array[int] = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
const DAYS_IN_EACH_MONTH : Array[int] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

## Date as a formatted string (DD-MM-YYYY).
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR | PROPERTY_USAGE_READ_ONLY)
var date_string: String = "":
	get(): return "%02d-%02d-%04d" % [day, month, year]

## Time as a formatted string (HH:MM:SS using a 24h clock).
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR | PROPERTY_USAGE_READ_ONLY)
var time_string: String = "":
	get(): return "%02d:%02d:%02d" % [hour, minute, second]

@export_group("Date")

## Day of the year.
@export
var day_of_the_year: int = 1:
	set(value):
		var day_from_zero = value - 1
		var current_year = year
		while day_from_zero < 0:
			current_year -= 1
			day_from_zero += days_in_year(current_year)
		while day_from_zero >= days_in_year(current_year):
			day_from_zero -= days_in_year(current_year)
			current_year += 1
		day_of_the_year = day_from_zero + 1
		year = current_year

### Day.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var day: int = 1:
	get():
		if is_leap_year(year):
			return day_of_the_year - CUMULATIVE_DAYS_IN_EACH_MONTH_IN_LEAP_YEAR[month - 1]
		else:
			return day_of_the_year - CUMULATIVE_DAYS_IN_EACH_MONTH[month - 1]
	set(value): day_of_the_year += (value - day)

## Month.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var month: int = 1:
	get():
		var zero_based_month = 0
		if is_leap_year(year):
			while zero_based_month < 12 and day_of_the_year > CUMULATIVE_DAYS_IN_EACH_MONTH_IN_LEAP_YEAR[zero_based_month]:
				zero_based_month += 1
		else:
			while zero_based_month < 12 and day_of_the_year > CUMULATIVE_DAYS_IN_EACH_MONTH[zero_based_month]:
				zero_based_month += 1
		return zero_based_month # we've overshot by one, so we get one based month again
	set(value):
		var current_day_of_the_month = day
		var zero_based_month = value - 1
		while zero_based_month < 0: zero_based_month += 12; year -= 1
		while zero_based_month > 11: year += 1; zero_based_month -= 12
		if is_leap_year(year):
			day_of_the_year = current_day_of_the_month + CUMULATIVE_DAYS_IN_EACH_MONTH_IN_LEAP_YEAR[zero_based_month]
		else:
			day_of_the_year = current_day_of_the_month + CUMULATIVE_DAYS_IN_EACH_MONTH[zero_based_month]

## Year.
@export
var year: int = 2025:
	set(value):
		if (day_of_the_year >= 60) and not (month == 2 and day == 60):
			if not is_leap_year(year) and is_leap_year(value):
				day_of_the_year += 1
			if is_leap_year(year) and not is_leap_year(value):
				day_of_the_year -= 1
		year = value

@export_group("Time")

## Time of the day in seconds.
@export
var time_in_seconds: float = 1:
	set(value):
		if value < 0.0 or value >= SECONDS_PER_DAY_F:
			var new_time := fposmod(value, SECONDS_PER_DAY_F)
			var delta_days := int((value - new_time) / SECONDS_PER_DAY_F)
			time_in_seconds = new_time
			if delta_days != 0:
				day_of_the_year += delta_days
		else:
			time_in_seconds = value
		#print(time_string)

## Hour.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var hour: int:
	get(): return int(time_in_seconds / SECONDS_PER_HOUR_F)
	set(value): time_in_seconds += (value - hour) * SECONDS_PER_HOUR_I

## Minute.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var minute: int:
	get(): return int(time_in_seconds / SECONDS_PER_MINUTE_F) % SECONDS_PER_MINUTE_I
	set(value): time_in_seconds += (value - minute) * SECONDS_PER_MINUTE_I

## Second.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var second: int = 1:
	get(): return int(time_in_seconds) % SECONDS_PER_MINUTE_I
	set(value): time_in_seconds += value - second

## Fraction of a second.
@export_custom(PROPERTY_HINT_NONE, "", PROPERTY_USAGE_EDITOR)
var second_fraction: float:
	get(): return time_in_seconds - int(time_in_seconds)
	set(value): time_in_seconds += value - second_fraction

## UTC offset in hours.
@export
var utc_offset_in_hours: float = 2.0

@export_group("Location")

## Latitude.
@export_custom(PROPERTY_HINT_NONE, "suffix:°") 
var latitude: float

## Longitude.
@export_custom(PROPERTY_HINT_NONE, "suffix:°") 
var longitude: float

## Altitude.
@export_custom(PROPERTY_HINT_NONE, "suffix:m")
var altitude: float

@export_group("Simulation")

## The speed at which time runs in this clock.
@export_range(0.0, 3600.0, 0.1)
var time_factor: float = 1.0

## Advance time while in the editor?
@export
var advance_time_in_editor: bool = false

## Sets the current date and time to the wall clock as configured on the local system.
@export_tool_button("Set date and time to wall clock")
var set_time_action = set_date_and_time_to_wall_clock

var cloud_sky := load(get_script().resource_path.get_base_dir() + "/clouds_sky.tres")

func _ready():
	cloud_sky.date_time_location = self

func _process(delta: float) -> void:
	if not advance_time_in_editor and Engine.is_editor_hint():
		return
	time_in_seconds += delta * time_factor

## Returns true if the provided year is a leap year.
func is_leap_year(year: int) -> bool:
	return (year % 4 == 0) and (year % 400 == 0 or year % 100 != 0)

## Returns the number of days in the provided year, taking into account leap years.
func days_in_year(year: int) -> int:
	return 366 if is_leap_year(year) else 365

## Sets the date and time to the wall clock as configured on the local system.
func set_date_and_time_to_wall_clock() -> void:
	var now := Time.get_datetime_dict_from_system()
	var timezone := Time.get_time_zone_from_system()
	year = now.year
	month = now.month
	day = now.day
	hour = now.hour
	minute = now.minute
	second = now.second
	utc_offset_in_hours = float(timezone.bias) / MINUTES_PER_HOUR_F

## Converts the local date and time to UTC and converts time in seconds to time in hours.
func convert_to_utc_date_and_time():
	var day_of_the_year_utc = day_of_the_year
	var year_utc = year
	var time_in_hours_utc = SECONDS_TO_HOURS_F * time_in_seconds - utc_offset_in_hours
	if time_in_hours_utc < 0.0:
		day_of_the_year_utc -= 1
		if day_of_the_year_utc < 0:
			year_utc -= 1
			day_of_the_year_utc = days_in_year(year_utc)
	if time_in_hours_utc >= 24.0:
		day_of_the_year_utc += 1
		if day_of_the_year_utc > days_in_year(year_utc):
			year_utc += 1
			day_of_the_year_utc = 1
	return {"time_in_hours": time_in_hours_utc, "day_of_the_year": day_of_the_year_utc, "year": year_utc}
