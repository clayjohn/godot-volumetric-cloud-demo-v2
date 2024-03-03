@tool
extends Texture2DRD

# This generates an LUT when it is loaded that is used throughout the sky system

var texture_size := Vector2i(256, 64)

var rd : RenderingDevice

var shader : RID
var pipeline : RID
var texture_rd : RID
var texture_set : RID

func _init():
	rd = RenderingServer.get_rendering_device()
	RenderingServer.call_on_render_thread(_initialize_texture)
	RenderingServer.call_on_render_thread.call_deferred(_initialize_compute_code)

func _notification(what):
	if what == NOTIFICATION_PREDELETE:
		if texture_rd:
			rd.free_rid(texture_rd)
		if shader:
			rd.free_rid(shader)

func _create_uniform_set(p_texture_rd : RID) -> RID:
	var uniform := RDUniform.new()
	uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	uniform.binding = 0
	uniform.add_id(p_texture_rd)
	return rd.uniform_set_create([uniform], shader, 0)
	

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

	# Create our texture
	texture_rd = rd.texture_create(tf, RDTextureView.new(), [])

func _initialize_compute_code():
	# Create our shader
	var shader_file = load(get_script().resource_path.get_base_dir() + "/transmittance-lut.glsl")
	var shader_spirv: RDShaderSPIRV = shader_file.get_spirv()
	shader = rd.shader_create_from_spirv(shader_spirv)
	if not shader.is_valid():
		print("you got an invalid shader buddy")
		return
	pipeline = rd.compute_pipeline_create(shader)

	# Now create our uniform set so we can use these textures in our shader
	texture_set = _create_uniform_set(texture_rd)

	texture_rd_rid = texture_rd

	var push_constant : PackedFloat32Array = PackedFloat32Array()
	push_constant.push_back(texture_size.x)
	push_constant.push_back(texture_size.y)
	push_constant.push_back(0.0)
	push_constant.push_back(0.0)

	# Run our compute shader.
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, pipeline)
	rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
	rd.compute_list_bind_uniform_set(compute_list, texture_set, 0)
	rd.compute_list_dispatch(compute_list, 32, 8, 1)
	rd.compute_list_end()
