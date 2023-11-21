shader_type sky;

uniform sampler2D hemisphere;

uniform sampler2D blend_from_texture : filter_linear, repeat_disable, source_color;
uniform sampler2D blend_to_texture : filter_linear, repeat_disable, source_color;
uniform float blend_amount : hint_range(0.0, 1.0, 0.01);
uniform vec4 ground_color : source_color = vec4(0.5);

void sky() {
	vec3 norm = EYEDIR;
	norm.y = max(0.0, norm.y);
	norm = normalize(norm);
	norm.xz /= 1.0 + abs(norm.y);

	// If facing down, clamp to horizon, then expand slightly.
	if (EYEDIR.y < 0.0) {
		norm.xz *= 1.05;
	}
	
	norm.xz = norm.xz * vec2(0.5, 0.5) + vec2(0.5, 0.5);
	vec3 blend_from = textureLod(blend_from_texture, norm.xz, 0.0).rgb;
	vec3 blend_to = textureLod(blend_to_texture, norm.xz, 0.0).rgb;
	COLOR = mix(blend_from, blend_to, blend_amount);
}
