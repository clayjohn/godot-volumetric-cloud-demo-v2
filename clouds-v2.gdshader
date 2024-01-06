shader_type sky;

uniform sampler2D hemisphere;

uniform sampler2D blend_from_texture : filter_linear, repeat_disable, source_color;
uniform sampler2D blend_to_texture : filter_linear, repeat_disable, source_color;
uniform sampler2D sky_texture : filter_linear, repeat_disable, source_color;
uniform float blend_amount : hint_range(0.0, 1.0, 0.01);
uniform vec4 ground_color : source_color = vec4(0.5);

vec2 oct_wrap(vec2 v) {
	vec2 signVal;
	signVal.x = v.x >= 0.0 ? 1.0 : -1.0;
	signVal.y = v.y >= 0.0 ? 1.0 : -1.0;
	return (1.0 - abs(v.yx)) * signVal;
}

vec2 vec3_to_oct(vec3 e) {
	e /= abs(e.x) + abs(e.y) + abs(e.z);
	e.xy = e.z >= 0.0 ? e.xy : oct_wrap( e.xy );

	vec2 n;
	n.y = e.y * 0.5 + 0.5;
	n.x = e.x * 0.5 + n.y;
	n.y = e.x * -0.5 + n.y;

	return n;
}

void sky() {
	vec3 norm = EYEDIR;
	norm.y = max(0.0, norm.y);
	norm = normalize(norm);

	norm.xz = vec3_to_oct(norm.xzy);
	vec3 blend_from = textureLod(blend_from_texture, norm.xz, 0.0).rgb;
	vec3 blend_to = textureLod(blend_to_texture, norm.xz, 0.0).rgb;
	COLOR = mix(blend_from, blend_to, blend_amount);
}
