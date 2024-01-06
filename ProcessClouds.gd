@tool
extends Node3D

class FrameData:
	# Internal properties
	var texture_size : int = 0
	var update_position : Vector2i = Vector2i(0, 0)

	# User-configurable properties
	var wind_direction : Vector2 = Vector2(1.0, 0.0)
	var wind_speed : float = 1.0
	var density : float = 0.05
	var cloud_coverage : float = 0.25
	var time_offset : float = 0.0
	var turbidity : float = 10.0
	var sun_disk_scale : float = 1.0
	var ground_color : Color = Color(1.0, 1.0, 1.0, 1.0)
	var exposure : float = 0.1
	var sun_visible : bool = true
	var moon_visible : bool = false

	# Properties updated by the light
	var LIGHT_DIRECTION : Vector3 = Vector3(0.0, 1.0, 0.0)
	var LIGHT_ENERGY : float = 1.0
	var LIGHT_COLOR : Color = Color(1.0, 1.0, 1.0, 1.0)
	
	func update_light_data(light : Light3D):
		LIGHT_DIRECTION = (light.transform.basis * Vector3(0.0, 0.0, 1.0)).normalized()
		LIGHT_ENERGY = light.light_energy
		LIGHT_COLOR = light.light_color.srgb_to_linear()

@export
var light_path : NodePath
@export
var wind_direction : Vector2 = Vector2(1.0, 0.0)
@export
var wind_speed : float = 1.0
@export
var density : float = 0.05
@export
var cloud_coverage : float = 0.25
@export
var time_offset : float = 0.0
@export
var turbidity : float = 10.0
@export
var sun_disk_scale : float = 1.0
@export
var ground_color : Color = Color(1.0, 1.0, 1.0, 1.0)
@export
var exposure : float = 0.1
@export
var sun_visible : bool = true
@export
var moon_visible : bool = false

var frame_data : FrameData = FrameData.new()
# Needs to be divisble by sqrt(frames_to_update)
var texture_size : int = 768
# needs to always be a power of two value
var frames_to_update : int = 64
var update_region_size : int = 96
var num_workgroups : int = 12

@export
var textures : Array = []
var texture_to_update : int = 0
var texture_to_blend_from : int = 1
var texture_to_blend_to : int = 2

var frame = 0

var can_run = false

func _ready():
	frame_data.texture_size = texture_size
	frame_data.update_light_data(get_node(light_path))
	# In case we're running stuff on the rendering thread
	# we need to do our initialization on that thread
	RenderingServer.call_on_render_thread(_initialize_compute_code.bind(texture_size))

	var frames_sqrt = sqrt(frames_to_update)
	update_region_size = texture_size / frames_sqrt
	num_workgroups = update_region_size / 8

	_update_per_frame_data()

	for i in range(frames_to_update * 2):
		update_sky()

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	if len(textures) < 3:
		return
		
	if not can_run:
		return
	
	update_sky()
	$Label.text = "Frame Time: " + str(1000.0/Engine.get_frames_per_second()) + " ms"
	
func update_sky():
	if frame >= frames_to_update:
		
		# Increase our next texture index
		texture_to_update = (texture_to_update + 1) % 3
		texture_to_blend_from = (texture_to_blend_from + 1) % 3
		texture_to_blend_to = (texture_to_blend_to + 1) % 3
		
		_update_per_frame_data() # Only call once per update otherwise quads get out of sync

		$WorldEnvironment.environment.sky.sky_material.set_shader_parameter("blend_from_texture", textures[texture_to_blend_from])
		$WorldEnvironment.environment.sky.sky_material.set_shader_parameter("blend_to_texture", textures[texture_to_blend_to])
		$WorldEnvironment.environment.sky.sky_material.set_shader_parameter("ground_color", frame_data.ground_color)
		$MeshInstance3D.material_override.albedo_texture = textures[texture_to_blend_to]
		frame = 0

	$WorldEnvironment.environment.sky.sky_material.set_shader_parameter("blend_amount", float(frame) / float(frames_to_update))

	RenderingServer.call_on_render_thread(_render_process.bind(texture_to_update))
	
	frame_data.update_position.x += update_region_size
	if frame_data.update_position.x >= texture_size:
		frame_data.update_position.x = 0
		frame_data.update_position.y += update_region_size
	if frame_data.update_position.y >= texture_size:
		frame_data.update_position = Vector2i(0, 0)
		
	frame += 1

###############################################################################
# Everything after this point is designed to run on our rendering thread

var rd : RenderingDevice

var shader : RID
var pipeline : RID

# We use 3 textures:
# - One to render into
# - One that contains the last frame rendered
# - One for the frame before that
var texture_rd : Array = [ RID(), RID(), RID() ]
var texture_set : Array = [ RID(), RID(), RID() ]
var noise_uniform_set : RID = RID()
var sampler

func _create_uniform_set(texture_rd : RID) -> RID:
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	uniform.binding = 0
	uniform.add_id(texture_rd)
	return rd.uniform_set_create([uniform], shader, 0)
	
func _create_noise_uniform_set() -> RID:
	var uniforms = []
	
	var sampler_state = RDSamplerState.new()
	sampler_state.repeat_u = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.repeat_v = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.repeat_w = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.mip_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	
	sampler = rd.sampler_create(sampler_state)
	
	var large_scale_noise = preload("res://perlworlnoise.tga")
	var LSN_rd = RenderingServer.texture_get_rd_texture(large_scale_noise.get_rid())
	
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 0
	uniform.add_id(sampler)
	uniform.add_id(LSN_rd)
	uniforms.push_back(uniform)
	
	var small_scale_noise = preload("res://worlnoise.bmp")
	var SSN_rd = RenderingServer.texture_get_rd_texture(small_scale_noise.get_rid())
	
	uniform = RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 1
	uniform.add_id(sampler)
	uniform.add_id(SSN_rd)
	uniforms.push_back(uniform)
	
	var weather_noise = preload("res://weather.bmp")
	var W_rd = RenderingServer.texture_get_rd_texture(weather_noise.get_rid())
	
	uniform = RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 2
	uniform.add_id(sampler)
	uniform.add_id(W_rd)
	uniforms.push_back(uniform)
	
	print(uniforms)
	
	return rd.uniform_set_create(uniforms, shader, 1)

func _initialize_compute_code(p_texture_size):
	can_run = true
	# As this becomes part of our normal frame rendering,
	# we use our main rendering device here.
	rd = RenderingServer.get_rendering_device()

	# Create our shader
	var shader_file = load("res://clouds.glsl")
	var shader_spirv: RDShaderSPIRV = shader_file.get_spirv()
	shader = rd.shader_create_from_spirv(shader_spirv)
	if not shader.is_valid():
		can_run = false
	pipeline = rd.compute_pipeline_create(shader)

	# Create our textures to manage our wave
	var tf : RDTextureFormat = RDTextureFormat.new()
	tf.format = RenderingDevice.DATA_FORMAT_R16G16B16A16_SFLOAT
	tf.texture_type = RenderingDevice.TEXTURE_TYPE_2D
	tf.width = p_texture_size
	tf.height = p_texture_size
	tf.depth = 1
	tf.array_layers = 1
	tf.mipmaps = 1
	tf.usage_bits = RenderingDevice.TEXTURE_USAGE_SAMPLING_BIT + RenderingDevice.TEXTURE_USAGE_COLOR_ATTACHMENT_BIT + RenderingDevice.TEXTURE_USAGE_STORAGE_BIT + RenderingDevice.TEXTURE_USAGE_CAN_UPDATE_BIT + RenderingDevice.TEXTURE_USAGE_CAN_COPY_TO_BIT

	noise_uniform_set = _create_noise_uniform_set()

	for i in range(3):
		# Create our texture
		texture_rd[i] = rd.texture_create(tf, RDTextureView.new(), [])

		# Make sure our textures are cleared
		rd.texture_clear(texture_rd[i], Color(float(i==0), float(i==1), float(i==2), 0), 0, 1, 0, 1)

		# Now create our uniform set so we can use these textures in our shader
		texture_set[i] = _create_uniform_set(texture_rd[i])
		textures.push_back(Texture2DRD.new())
		#textures[i].texture_rd_rid = texture_rd[i]

func _update_per_frame_data():
	frame_data.update_light_data(get_node(light_path))
	frame_data.wind_direction = wind_direction
	frame_data.wind_speed = wind_speed
	frame_data.density = density
	frame_data.cloud_coverage = cloud_coverage
	frame_data.time_offset = time_offset
	frame_data.turbidity = turbidity
	frame_data.sun_disk_scale = sun_disk_scale
	frame_data.ground_color = ground_color
	frame_data.exposure = exposure
	
func _fill_push_constant():
	
	var push_constant : PackedFloat32Array = PackedFloat32Array()
#	The order of properties here must match those in clouds.glsl, including padding
	push_constant.push_back(frame_data.texture_size)
	push_constant.push_back(frame_data.texture_size)
	push_constant.push_back(frame_data.update_position.x)
	push_constant.push_back(frame_data.update_position.y)

	push_constant.push_back(frame_data.wind_direction.x)
	push_constant.push_back(frame_data.wind_direction.y)
	push_constant.push_back(frame_data.wind_speed)
	push_constant.push_back(frame_data.density)
	
	push_constant.push_back(frame_data.cloud_coverage)
	push_constant.push_back(frame_data.time_offset)
	push_constant.push_back(frame_data.turbidity)
	push_constant.push_back(frame_data.sun_disk_scale)
	
	push_constant.push_back(frame_data.ground_color.r)
	push_constant.push_back(frame_data.ground_color.g)
	push_constant.push_back(frame_data.ground_color.b)
	push_constant.push_back(frame_data.ground_color.a)
	
	push_constant.push_back(frame_data.LIGHT_DIRECTION.x)
	push_constant.push_back(frame_data.LIGHT_DIRECTION.y)
	push_constant.push_back(frame_data.LIGHT_DIRECTION.z)
	push_constant.push_back(frame_data.LIGHT_ENERGY)
	
	push_constant.push_back(frame_data.LIGHT_COLOR.r)
	push_constant.push_back(frame_data.LIGHT_COLOR.g)
	push_constant.push_back(frame_data.LIGHT_COLOR.b)
	push_constant.push_back(0.0) # will be flags
	
	push_constant.push_back(frame_data.exposure)
	push_constant.push_back(Time.get_ticks_msec()/1000.0)
	push_constant.push_back(0.0)
	push_constant.push_back(0.0)

	return push_constant

func _render_process(p_texture_to_update):
	
	#if not textures[p_texture_to_update].is_valid():
	textures[p_texture_to_update].texture_rd_rid = texture_rd[p_texture_to_update]

	var push_constant = _fill_push_constant()

	# Run our compute shader.
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, pipeline)
	rd.compute_list_bind_uniform_set(compute_list, noise_uniform_set, 1)
	rd.compute_list_bind_uniform_set(compute_list, texture_set[p_texture_to_update], 0)

	rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
	rd.compute_list_dispatch(compute_list, num_workgroups, num_workgroups, 1)
	rd.compute_list_end()

func _free_compute_resources():
	# Note that our sets and pipeline are cleaned up automatically as they are dependencies :P
	for i in range(3):
		if texture_rd[i]:
			rd.free_rid(texture_rd[i])

	if shader:
		rd.free_rid(shader)
