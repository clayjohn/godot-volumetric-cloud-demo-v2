@tool
extends DirectionalLight3D

var atmosphere_tex := load(get_script().resource_path.get_base_dir() + "/sky_lut.tres")
var cloud_sky := load(get_script().resource_path.get_base_dir() + "/clouds_sky.tres")

func _ready():
	call_deferred("delayed_init")
	
# Needs to be called after cloud_sky is fully initialized
func delayed_init():
	cloud_sky.sun = self
	cloud_sky.request_full_sky_init()
	
func _notification(what):
	if what == NOTIFICATION_TRANSFORM_CHANGED:
		atmosphere_tex.request_update()
