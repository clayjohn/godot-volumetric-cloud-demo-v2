shader_type sky;
render_mode use_debanding;

uniform sampler2D blend_from_texture : filter_linear, repeat_disable, source_color;
uniform sampler2D blend_to_texture : filter_linear, repeat_disable, source_color;
uniform sampler2D source_panorama : filter_linear, hint_default_black, repeat_disable;
uniform sampler2D source_transmittance : filter_linear, hint_default_black, repeat_disable;

uniform float blend_amount : hint_range(0.0, 1.0, 0.01);
uniform float sun_disk_scale = 1.0;

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

vec3 getValFromSkyLUT(vec3 rayDir) {
	vec2 uv;
	float phi = atan(rayDir.z, rayDir.x);
    float theta = asin(rayDir.y);
	uv.x = (phi / PI * 0.5 + 0.5);
    // Undo the non-linear transformation from the sky-view LUT
    uv.y = sqrt(abs(theta) / (PI * 0.5)) * sign(theta) * 0.5 + 0.5;
    return texture(source_panorama, uv).rgb/50.0;
}

// https://www.shadertoy.com/view/slSXRW
vec3 sunWithBloom(vec3 rayDir, vec3 sunDir) {
    float sunSolidAngle = sun_disk_scale * 0.53*PI/180.0;
    float minSunCosTheta = cos(sunSolidAngle);

    float cosTheta = dot(rayDir, sunDir);
    if (cosTheta >= minSunCosTheta) return vec3(1.0);
    
    float offset = minSunCosTheta - cosTheta;
    float gaussianBloom = exp(-offset*50000.0)*0.5;
    float invBloom = 1.0/(0.02 + offset*300.0)*0.01;
    return vec3(gaussianBloom+invBloom);
}

float rayIntersectSphere(vec3 ro, vec3 rd, float rad) {
    float b = dot(ro, rd);
    float c = dot(ro, ro) - rad*rad;
    if (c > 0.0f && b > 0.0) return -1.0;
    float discr = b*b - c;
    if (discr < 0.0) return -1.0;
    // Special case: inside sphere, use far discriminant
    if (discr > b*b) return (-b + sqrt(discr));
    return -b - sqrt(discr);
}

const float groundRadiusMM = 6.360;
const float atmosphereRadiusMM = 6.460;
const vec3 viewPos = vec3(0.0, groundRadiusMM + 0.0002, 0.0);
const vec2 tLUTRes = vec2(256.0, 64.0);

vec3 getValFromTLUT(sampler2D tex, vec2 bufferRes, vec3 pos, vec3 sunDir) {
    float height = length(pos);
    vec3 up = pos / height;
	float sunCosZenithAngle = dot(up,sunDir);
    vec2 uv = vec2(tLUTRes.x*clamp(0.5 + 0.5*sunCosZenithAngle, 0.0, 1.0),
                   tLUTRes.y*max(0.0, min(1.0, (height - groundRadiusMM)/(atmosphereRadiusMM - groundRadiusMM))));
    uv /= bufferRes;
    return texture(tex, uv).rgb;
}

vec3 get_atmo(vec3 dir) {
	vec3 col = getValFromSkyLUT(dir);
	vec3 sunLum = sunWithBloom(dir, LIGHT0_DIRECTION);
    // Use smoothstep to limit the effect, so it drops off to actual zero.
    sunLum = smoothstep(0.002, 1.0, sunLum);
    if (length(sunLum) > 0.0) {
        if (rayIntersectSphere(viewPos, dir, groundRadiusMM) >= 0.0) {
            sunLum *= 0.0;
        } else {
            // If the sun value is applied to this pixel, we need to calculate the transmittance to obscure it.
        	sunLum *= getValFromTLUT(source_transmittance, tLUTRes, viewPos, LIGHT0_DIRECTION);
		}
    }
    col += sunLum;
	return col;
}

void sky() {
	vec3 norm = EYEDIR;
	norm.y = max(0.0, norm.y);
	norm = normalize(norm);

	norm.xz = vec3_to_oct(norm.xzy);
	vec4 blend_from = textureLod(blend_from_texture, norm.xz, 0.0);
	vec4 blend_to = textureLod(blend_to_texture, norm.xz, 0.0);
	vec4 clouds = mix(blend_from, blend_to, blend_amount);
	vec3 background = get_atmo(EYEDIR);
	COLOR = background * (1.0-clouds.a) + clouds.xyz;
	COLOR = mix(clamp(COLOR, vec3(0.0), vec3(100.0)), clamp(background, vec3(0.0), vec3(100.0)), smoothstep(0.6, 1.0, 1.0 - EYEDIR.y));
}
