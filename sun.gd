@tool
extends DirectionalLight3D

@export
var atmosphere_tex : Texture2DRD:
	set(tex):
		atmosphere_tex = tex
		atmosphere_tex.sun = self

@export
var world_environment : WorldEnvironment

func _ready():
	world_environment.environment.sky.sun = self
	world_environment.environment.sky.initialize_sky()
	atmosphere_tex.sun = self
	
func _notification(what):
	if what == NOTIFICATION_TRANSFORM_CHANGED:
		atmosphere_tex.sun = self
		atmosphere_tex.update_lut()
