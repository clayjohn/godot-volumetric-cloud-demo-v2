@tool
extends Texture2DRD

var texture_size := Vector2i(200, 100)
var light_direction := Vector3(0.0, -1.0, 0.0)
var needs_update = true
var initialized = false

var back_texture : Array = [ Texture2DRD.new(), Texture2DRD.new() ]
var needs_full_update = true

var rd : RenderingDevice
var shader : RID
var pipeline : RID
# We need three copies so we can synchronize with the clouds
var texture_rd : Array = [ RID(), RID(), RID() ]
var texture_set : Array = [ RID(), RID(), RID() ]
var current_texture = 0
var lut_set : RID

var transmittance_tex := load(get_script().resource_path.get_base_dir() + "/transmittance_lut.tres")

func _init():
	rd = RenderingServer.get_rendering_device()
	RenderingServer.call_on_render_thread(_initialize_texture)
	RenderingServer.call_on_render_thread.call_deferred(_initialize_compute_code)

func request_update():
	needs_update = true

# This should be called once to update the entire sky.
func update_lut(sun_direction : Vector3):
	light_direction = sun_direction
	if not initialized:
		print("Attempting to update uninitialized sky lut")
		return
	RenderingServer.call_on_render_thread(render_lut)
	if needs_full_update:
		RenderingServer.call_on_render_thread(render_lut)
		RenderingServer.call_on_render_thread(render_lut)
		needs_full_update = false

func _create_uniform_set(p_texture_rd : RID) -> RID:
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	uniform.binding = 0
	uniform.add_id(p_texture_rd)
	return rd.uniform_set_create([uniform], shader, 0)
	
func _create_uniform_set_for_sampling() -> RID:
	var sampler_state = RDSamplerState.new()
	sampler_state.repeat_u = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.repeat_v = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.repeat_w = RenderingDevice.SAMPLER_REPEAT_MODE_CLAMP_TO_EDGE
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.mip_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	
	var sampler = rd.sampler_create(sampler_state)
	
	var uniforms = []
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
	uniform.binding = 0
	uniform.add_id(sampler)
	uniform.add_id(transmittance_tex.texture_rd)
	uniforms.push_back(uniform)
	
	return rd.uniform_set_create(uniforms, shader, 1)

func _initialize_texture():
	var tf : RDTextureFormat = RDTextureFormat.new()
	tf.format = RenderingDevice.DATA_FORMAT_R16G16B16A16_SFLOAT
	tf.texture_type = RenderingDevice.TEXTURE_TYPE_2D
	tf.width = texture_size.x
	tf.height = texture_size.y
	tf.depth = 1
	tf.array_layers = 1
	tf.mipmaps = 1
	tf.usage_bits = RenderingDevice.TEXTURE_USAGE_SAMPLING_BIT + RenderingDevice.TEXTURE_USAGE_COLOR_ATTACHMENT_BIT + RenderingDevice.TEXTURE_USAGE_STORAGE_BIT + RenderingDevice.TEXTURE_USAGE_CAN_UPDATE_BIT + RenderingDevice.TEXTURE_USAGE_CAN_COPY_TO_BIT
	if Engine.is_editor_hint():
		tf.usage_bits += RenderingDevice.TEXTURE_USAGE_CAN_COPY_FROM_BIT

	# Create our textures
	texture_rd[0] = rd.texture_create(tf, RDTextureView.new(), [])
	texture_rd[1] = rd.texture_create(tf, RDTextureView.new(), [])
	texture_rd[2] = rd.texture_create(tf, RDTextureView.new(), [])

func _initialize_compute_code():
	# Create our shader
	var shader_file = load(get_script().resource_path.get_base_dir() + "/sky-lut.glsl")
	var shader_spirv: RDShaderSPIRV = shader_file.get_spirv()
	shader = rd.shader_create_from_spirv(shader_spirv)
	if not shader.is_valid():
		print("you got an invalid shader buddy")
		return
	pipeline = rd.compute_pipeline_create(shader)

	# Now create our uniform set so we can use these textures in our shader
	texture_set[0] = _create_uniform_set(texture_rd[0])
	texture_set[1] = _create_uniform_set(texture_rd[1])
	texture_set[2] = _create_uniform_set(texture_rd[2])
	lut_set = _create_uniform_set_for_sampling()

	texture_rd_rid = texture_rd[0]
	back_texture[0].texture_rd_rid = texture_rd[1]
	back_texture[1].texture_rd_rid = texture_rd[2]
	
	initialized = true

func render_lut():
	var push_constant : PackedFloat32Array = PackedFloat32Array()
	push_constant.push_back(texture_size.x)
	push_constant.push_back(texture_size.y)
	push_constant.push_back(0.0)
	push_constant.push_back(0.0)
	
	push_constant.push_back(light_direction.x)
	push_constant.push_back(light_direction.y)
	push_constant.push_back(light_direction.z)
	push_constant.push_back(0.0)

	# Run our compute shader.
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, pipeline)
	rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
	rd.compute_list_bind_uniform_set(compute_list, texture_set[current_texture], 0)
	rd.compute_list_bind_uniform_set(compute_list, lut_set, 1)
	rd.compute_list_dispatch(compute_list, 25, 13, 1)
	rd.compute_list_end()

	texture_rd_rid = texture_rd[current_texture]
	current_texture = (current_texture + 1) % 3
	back_texture[0].texture_rd_rid = texture_rd[current_texture]
	back_texture[1].texture_rd_rid = texture_rd[(current_texture + 1) % 3]
	
	needs_update = false
