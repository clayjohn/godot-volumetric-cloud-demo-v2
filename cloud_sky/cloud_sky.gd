@tool
extends Sky

# User configurable properties
@export_group("Cloud Settings")

# Sets the wind direction in degrees, where 0 degrees is a wind coming from the north (+X axis),
# 90 degrees is east, 180 is south and 270 (or rather -90) is west.
@export_custom(PROPERTY_HINT_RANGE, "-180,180,0.1,radians_as_degrees")
var wind_direction: float = 0.0

# Sets the wind speed in m/s.
@export_custom(PROPERTY_HINT_RANGE, "0,120,0.1,or_greater,or_less,suffix:m/s")
var wind_speed: float = 1.0
# TODO This needs calibration if we want the m/s to make some sense. We can't really specify the
#      altitude of the clouds right now so it's more about getting things in a believable range
#      where 120 m/s is quite a strong hurricane.

@export
var density : float = 0.05
@export
var cloud_coverage : float = 0.25
@export
var time_offset : float = 0.0

@export_group("Sky Settings")
@export
var sun_disk_scale : float = 1.0:
	set(scale):
		sun_disk_scale = scale
		sky_material.set_shader_parameter("sun_disk_scale", sun_disk_scale)
@export
var ground_color : Color = Color(1.0, 1.0, 1.0, 1.0)

@export_group("Performance Settings")
@export_enum("Very Fast(4):4", "Fast(16):16", "Default(64):64", "Performance(256):256")
var frames_to_update : int = 64:
	set(value):
		frames_to_update = value
		cleanup()
		update_performance()
		request_full_sky_init()

@export_range(32.0, 8192.0, 32.0)
var texture_size : int = 768: # Needs to be divisible by sqrt(frames_to_update)
	set(value):
		texture_size = value
		cleanup()
		update_performance()
		request_full_sky_init()

var sun : DirectionalLight3D

# Everything in the compute shader must be cached here so that it only updates
# after swapping to a new texture.
class FrameData:
	# User-configurable properties
	var wind_direction : Vector2 = Vector2(1.0, 0.0)
	var wind_speed : float = 1.0
	var density : float = 0.05
	var cloud_coverage : float = 0.25
	var time_offset : float = 0.0
	var ground_color : Color = Color(1.0, 1.0, 1.0, 1.0)
	
	# Properties we calculate based on the user-configurable ones
	var _time : float = 0.0
	var _cloud_pos : Vector2 = Vector2(0.0, 0.0)
	var _detailed_pos : Vector2 = Vector2(0.0, 0.0)
	var _weather_pos : Vector2 = Vector2(0.0, 0.0)

	# Properties updated by the light
	var LIGHT_DIRECTION : Vector3 = Vector3(0.0, -1.0, 0.0)
	var LIGHT_ENERGY : float = 1.0
	var LIGHT_COLOR : Color = Color(1.0, 1.0, 1.0, 1.0)
	
	func update_light_data(light : Light3D):
		LIGHT_DIRECTION = (light.transform.basis * Vector3(0.0, 0.0, 1.0)).normalized()
		LIGHT_ENERGY = light.light_energy
		LIGHT_COLOR = light.light_color.srgb_to_linear()

var frame_data : FrameData = FrameData.new()
var update_position : Vector2i = Vector2i(0, 0)
var update_region_size : int = 96 # texture_size / sqrt(frames_to_update)
var num_workgroups : int = 12 # update_region_size / 8

var textures : Array = []
var texture_to_update : int = 0
var texture_to_blend_from : int = 1
var texture_to_blend_to : int = 2

var sky_lut := load(get_script().resource_path.get_base_dir() + "/sky_lut.tres")
var transmittance_tex := load(get_script().resource_path.get_base_dir() + "/transmittance_lut.tres")

var frame = 0

var can_run = false
var needs_full_sky_init = true

func _init():
	call_deferred("delayed_init")
	
# Workaround due to the fact that exports are set after _init() is called
func delayed_init():
	update_performance()

	# This calls "update_sky" at the beginning of the render loop automatically.
	RenderingServer.connect("frame_pre_draw", update_sky)

func update_performance():
	var frames_sqrt : int = int(sqrt(frames_to_update))
	update_region_size = texture_size / frames_sqrt
	if texture_size % frames_sqrt !=0:
		texture_size = update_region_size * frames_sqrt
		print("texture_size is not a multiple of sqrt(frames_to_update), changing to: ", texture_size)
	num_workgroups = (update_region_size + 7) / 8

	sky_material.set_shader_parameter("sun_disk_scale", sun_disk_scale)
	RenderingServer.call_on_render_thread.call_deferred(_initialize_compute_code.bind(texture_size))

func request_full_sky_init():
	needs_full_sky_init = true

# Initialize and update the sky so it is visible in the first frame.
func initialize_sky():
	_update_per_frame_data()
	for i in range(frames_to_update * 2):
		update_sky()

func update_sky():
	if not can_run:
		return

	if needs_full_sky_init:
		needs_full_sky_init = false
		initialize_sky()

	if frame >= frames_to_update:
		# Increase our next texture index
		texture_to_update = (texture_to_update + 1) % 3
		texture_to_blend_from = (texture_to_blend_from + 1) % 3
		texture_to_blend_to = (texture_to_blend_to + 1) % 3
		_update_per_frame_data() # Only call once per update otherwise quads get out of sync

		sky_material.set_shader_parameter("blend_from_texture", textures[texture_to_blend_from])
		sky_material.set_shader_parameter("blend_to_texture", textures[texture_to_blend_to])
		
		sky_material.set_shader_parameter("sky_blend_from_texture", sky_lut.back_texture[0])
		sky_material.set_shader_parameter("sky_blend_to_texture", sky_lut.back_texture[1])

		frame = 0

	sky_material.set_shader_parameter("blend_amount", float(frame) / float(frames_to_update))

	RenderingServer.call_on_render_thread(_render_process.bind(texture_to_update))
	
	update_position.x += update_region_size
	if update_position.x >= texture_size:
		update_position.x = 0
		update_position.y += update_region_size
	if update_position.y >= texture_size:
		update_position = Vector2i(0, 0)
		
	frame += 1

func _update_per_frame_data():
	if sun:
		frame_data.update_light_data(sun)
	frame_data.wind_direction = Vector2.from_angle(wind_direction)
	frame_data.wind_speed = wind_speed
	frame_data.density = density
	frame_data.cloud_coverage = cloud_coverage
	frame_data.time_offset = time_offset
	frame_data.ground_color = ground_color

	# Before we can push those constants, there's a bit of math we need to do
	var time = Time.get_ticks_msec() / 1000.0
	var delta = time - frame_data._time
	var delta2 = delta * 0.001 + 0.005 * frame_data.time_offset;
	var wind_direction_normalized = frame_data.wind_direction.normalized()

	# Which involves keeping some state and integrating positions over the latest time step
	frame_data._time = time
	frame_data._detailed_pos += delta * wind_direction_normalized
	frame_data._cloud_pos += delta * wind_direction_normalized * frame_data.wind_speed
	frame_data._weather_pos += delta2 * wind_direction_normalized * frame_data.wind_speed

	sky_lut.update_lut(frame_data.LIGHT_DIRECTION)

func _validate_property(property):
	if property.name == "sky_material" or property.name == "process_mode":
		property.usage &= ~PROPERTY_USAGE_EDITOR

func _notification(what):
	if what == NOTIFICATION_PREDELETE:
		cleanup()
		
func cleanup():
	can_run = false
	frame = 0
	for i in range(3):
		if texture_rd[i]:
			rd.free_rid(texture_rd[i])
	if shader_rd:
		rd.free_rid(shader_rd)
	if noise_sampler:
		rd.free_rid(noise_sampler)
	if sky_sampler:
		rd.free_rid(sky_sampler)


###############################################################################
# Everything after this point is designed to run on our rendering thread

var rd : RenderingDevice

var shader_rd : RID
var pipeline : RID

# We use 3 textures:
# - One to render into
# - One that contains the last frame rendered
# - One for the frame before that
var texture_rd : Array = [ RID(), RID(), RID() ]
var texture_set : Array = [ RID(), RID(), RID() ]
var noise_uniform_set : RID = RID()
var sky_uniform_set : Array = [ RID(), RID(), RID() ]
var noise_sampler
var sky_sampler

func _render_process(p_texture_to_update):
	textures[p_texture_to_update].texture_rd_rid = texture_rd[p_texture_to_update]

	var push_constant = _fill_push_constant()

	# Run our compute shader.
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, pipeline)
	rd.compute_list_bind_uniform_set(compute_list, sky_uniform_set[(sky_lut.current_texture + 2) % 3], 2)
	rd.compute_list_bind_uniform_set(compute_list, noise_uniform_set, 1)
	rd.compute_list_bind_uniform_set(compute_list, texture_set[p_texture_to_update], 0)

	rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
	rd.compute_list_dispatch(compute_list, num_workgroups, num_workgroups, 1)
	rd.compute_list_end()

	
func _fill_push_constant():
	var push_constant : PackedFloat32Array = PackedFloat32Array()
	# The order of properties here must match those in clouds.glsl, including padding
	push_constant.push_back(texture_size)
	push_constant.push_back(texture_size)
	push_constant.push_back(update_position.x)
	push_constant.push_back(update_position.y)
	
	push_constant.push_back(frame_data._cloud_pos.x)
	push_constant.push_back(frame_data._cloud_pos.y)
	push_constant.push_back(frame_data._detailed_pos.x)
	push_constant.push_back(frame_data._detailed_pos.y)

	push_constant.push_back(frame_data._weather_pos.x)
	push_constant.push_back(frame_data._weather_pos.y)
	push_constant.push_back(0.0) # vec2 pad1
	push_constant.push_back(0.0) #
	
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
	push_constant.push_back(frame_data._time)
	
	push_constant.push_back(0.0) # float pad2
	push_constant.push_back(frame_data.density)
	push_constant.push_back(frame_data.cloud_coverage)
	push_constant.push_back(frame_data.time_offset)

	return push_constant

func _create_uniform_set(p_texture_rd : RID) -> RID:
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	uniform.binding = 0
	uniform.add_id(p_texture_rd)
	return rd.uniform_set_create([uniform], shader_rd, 0)

func _create_noise_uniform_set() -> RID:
	var uniforms = []
	
	var sampler_state = RDSamplerState.new()
	sampler_state.repeat_u = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.repeat_v = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.repeat_w = RenderingDevice.SAMPLER_REPEAT_MODE_REPEAT
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.mip_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	
	noise_sampler = rd.sampler_create(sampler_state)
	
	var large_scale_noise = preload("perlworlnoise.tga")
	var LSN_rd = RenderingServer.texture_get_rd_texture(large_scale_noise.get_rid())
	
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 0
	uniform.add_id(noise_sampler)
	uniform.add_id(LSN_rd)
	uniforms.push_back(uniform)
	
	var small_scale_noise = preload("worlnoise.bmp")
	var SSN_rd = RenderingServer.texture_get_rd_texture(small_scale_noise.get_rid())
	
	uniform = RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 1
	uniform.add_id(noise_sampler)
	uniform.add_id(SSN_rd)
	uniforms.push_back(uniform)
	
	var weather_noise = preload("weather.bmp")
	var W_rd = RenderingServer.texture_get_rd_texture(weather_noise.get_rid())
	
	uniform = RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 2
	uniform.add_id(noise_sampler)
	uniform.add_id(W_rd)
	uniforms.push_back(uniform)

	return rd.uniform_set_create(uniforms, shader_rd, 1)
	
func _create_sky_uniform_set(tex_id : int) -> RID:
	var uniforms = []
		
	var uniform = RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 0
	uniform.add_id(sky_sampler)
	uniform.add_id(sky_lut.texture_rd[tex_id])
	uniforms.push_back(uniform)
	
	return rd.uniform_set_create(uniforms, shader_rd, 2)

func _initialize_compute_code(p_texture_size):
	rd = RenderingServer.get_rendering_device()

	# Create our shader
	var shader_file = load(get_script().resource_path.get_base_dir() + "/clouds.glsl")
	var shader_spirv: RDShaderSPIRV = shader_file.get_spirv()
	shader_rd = rd.shader_create_from_spirv(shader_spirv)
	if not shader_rd.is_valid():
		can_run = false
		return
	pipeline = rd.compute_pipeline_create(shader_rd)

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
	if Engine.is_editor_hint():
		tf.usage_bits += RenderingDevice.TEXTURE_USAGE_CAN_COPY_FROM_BIT
	noise_uniform_set = _create_noise_uniform_set()

	var sampler_state = RDSamplerState.new()
	sampler_state = RDSamplerState.new()
	sampler_state.repeat_u = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.repeat_v = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.repeat_w = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.mip_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	
	sky_sampler = rd.sampler_create(sampler_state)

	sky_uniform_set[0] = _create_sky_uniform_set(0)
	sky_uniform_set[1] = _create_sky_uniform_set(1)
	sky_uniform_set[2] = _create_sky_uniform_set(2)
	textures.clear()

	for i in range(3):
		# Create our texture
		texture_rd[i] = rd.texture_create(tf, RDTextureView.new(), [])

		# Make sure our textures are cleared
		rd.texture_clear(texture_rd[i], Color(float(i==0), float(i==1), float(i==2), 0), 0, 1, 0, 1)

		# Now create our uniform set so we can use these textures in our shader
		texture_set[i] = _create_uniform_set(texture_rd[i])
		textures.push_back(Texture2DRD.new())

	can_run = true
