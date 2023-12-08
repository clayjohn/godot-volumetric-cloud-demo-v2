@tool
extends DirectionalLight3D

@export
var world_environment : WorldEnvironment

func _ready():
	world_environment.environment.sky.sun = self
	world_environment.environment.sky.initialize_sky()
