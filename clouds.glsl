#[compute]
#version 450

// Invocations in the (x, y, z) dimension
layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

// Our textures
layout(rgba16f, set = 0, binding = 0) uniform restrict writeonly image2D current_image;

layout(set = 1, binding = 0) uniform sampler3D large_scale_noise;
layout(set = 1, binding = 1) uniform sampler3D small_scale_noise;
layout(set = 1, binding = 2) uniform sampler2D weather_noise;
layout(set = 1, binding = 3) uniform sampler2D sky_lut;
layout(set = 1, binding = 4) uniform sampler2D transmittance_lut;

/*
TODOs ======================================
1. Add sun and moon
2. Optimize by using the old frame texture
3. (DONE)~Better noise to hide artifacts~
4. (DONE)~Use octahedral mapping (either to get the full sky, or just increase precision)~
5. Implement heirarchical marching.
6. (DONE)~Convert everything into a sky resource, or a sky material~
7. Write Atmosphere to a 64x64 cubemap and sample it here.
*/


#define SUN_ENABLED 1
#define MOON_ENABLED 2

// Our push constant.
// Push constants have a max size of 128 bytes (32 floats).
layout(push_constant, std430) uniform Params {
	vec2 texture_size;
	vec2 update_position;

	vec2 wind_direction;
	float wind_speed;
	float density;

	float cloud_coverage;
	float time_offset;
	float turbidity;
	float sun_disk_scale;

	vec4 ground_color;

	vec3 LIGHT_DIRECTION;
	float LIGHT_ENERGY;

	vec3 LIGHT_COLOR;
	uint flags;

	float exposure;
	float time;
	float custom_scale;
	float pad2;
} params;

const float rayleigh = 2.0;
const vec4 rayleigh_color = vec4(0.055, 0.14, 0.2957, 1.0);
const float mie = 0.005;
const float mie_eccentricity = 0.8;
const vec4 mie_color = vec4(0.3547, 0.5542, 0.8276, 1.0);

// Approximately earth sizes
const float g_radius = 6000000.0; //ground radius
const float sky_b_radius = 6001500.0;//bottom of cloud layer
const float sky_t_radius = 6004000.0;//top of cloud layer

const vec3 UP = vec3( 0.0, 1.0, 0.0 );

// Sun constants
const float SOL_SIZE = 0.00872663806;
const float SUN_ENERGY = 1000.0;

// optical length at zenith for molecules
const float rayleigh_zenith_size = 8.4e3;
const float mie_zenith_size = 1.25e3;

const float PI = 3.141592;

const float groundRadiusMM = 6.360;
const float atmosphereRadiusMM = 6.460;
const vec3 viewPos = vec3(0.0, groundRadiusMM + 0.0002, 0.0);
const vec2 tLUTRes = vec2(256.0, 64.0);

float safeacos(const float x) {
    return acos(clamp(x, -1.0, 1.0));
}

vec3 getValFromSkyLUT(vec3 rayDir) {
    float height = length(viewPos);
    vec3 up = viewPos / height;
    
    float horizonAngle = safeacos(sqrt(height * height - groundRadiusMM * groundRadiusMM) / height);
    float altitudeAngle = horizonAngle - acos(dot(rayDir, up)); // Between -PI/2 and PI/2
    float azimuthAngle; // Between 0 and 2*PI
    if (abs(altitudeAngle) > (0.5*PI - .0001)) {
        // Looking nearly straight up or down.
        azimuthAngle = 0.0;
    } else {
        vec3 right = cross(up, params.LIGHT_DIRECTION);
        vec3 forward = cross(up, right);
        
        vec3 projectedDir = normalize(rayDir - up*(dot(rayDir, up)));
        float sinTheta = dot(projectedDir, right);
        float cosTheta = dot(projectedDir, forward);
        azimuthAngle = atan(sinTheta, cosTheta) + PI;
    }
    
    // Non-linear mapping of altitude angle. See Section 5.3 of the paper.
    float v = 0.5 + 0.5*sign(altitudeAngle)*sqrt(abs(altitudeAngle)*2.0/PI);
    vec2 uv = vec2(azimuthAngle / (2.0*PI), v);
    
    return texture(sky_lut, uv).rgb;
}

vec3 sunWithBloom(vec3 rayDir, vec3 sunDir) {
    const float sunSolidAngle = 0.53*PI/180.0;
    const float minSunCosTheta = cos(sunSolidAngle);

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

vec3 getValFromTLUT(sampler2D tex, vec2 bufferRes, vec3 pos, vec3 sunDir) {
    float height = length(pos);
    vec3 up = pos / height;
	float sunCosZenithAngle = dot(up,sunDir);
    vec2 uv = vec2(tLUTRes.x*clamp(0.5 + 0.5*sunCosZenithAngle, 0.0, 1.0),
                   tLUTRes.y*max(0.0, min(1.0, (height - groundRadiusMM)/(atmosphereRadiusMM - groundRadiusMM))));
    uv /= bufferRes;
    return texture(tex, uv).rgb;
}

// From: https://www.shadertoy.com/view/4sfGzS credit to iq
float hash(vec3 p) {
	p  = fract( p * 0.3183099 + 0.1 );
	p *= 17.0;
	return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}

// Utility function that maps a value from one range to another. 
float remap(float originalValue,  float originalMin,  float originalMax,  float newMin,  float newMax) {
	return newMin + (((originalValue - originalMin) / (originalMax - originalMin)) * (newMax - newMin));
}

// Phase function
float henyey_greenstein(float cos_theta, float g) {
	const float k = 0.0795774715459;
	return k * (1.0 - g * g) / (pow(1.0 + g * g - 2.0 * g * cos_theta, 1.5));
}

vec3 get_sun(vec3 eye_dir) {
	vec3 sunLum = sunWithBloom(eye_dir, params.LIGHT_DIRECTION);
    // Use smoothstep to limit the effect, so it drops off to actual zero.
    sunLum = smoothstep(0.002, 1.0, sunLum);
    if (length(sunLum) > 0.0) {
        if (rayIntersectSphere(viewPos, eye_dir, groundRadiusMM) >= 0.0) {
            sunLum *= 0.0;
        } else {
            // If the sun value is applied to this pixel, we need to calculate the transmittance to obscure it.
        	sunLum *= getValFromTLUT(transmittance_lut, tLUTRes, viewPos, params.LIGHT_DIRECTION);
		}
    }
	return sunLum * 10.0;
}

// Simple Analytic sky. In a real project you should use a texture
vec3 atmosphere(vec3 eye_dir) {
	vec3 COLOR = getValFromSkyLUT(eye_dir);
	vec3 sunLum = get_sun(eye_dir);
	COLOR *= 10.0;
    COLOR += sunLum;
	return COLOR;
}

float GetHeightFractionForPoint(float inPosition) { 
	float height_fraction = (inPosition -  sky_b_radius) / (sky_t_radius - sky_b_radius); 
	return clamp(height_fraction, 0.0, 1.0);
}

vec4 mixGradients(float cloudType){
	const vec4 STRATUS_GRADIENT = vec4(0.02f, 0.05f, 0.09f, 0.11f);
	const vec4 STRATOCUMULUS_GRADIENT = vec4(0.02f, 0.2f, 0.48f, 0.625f);
	const vec4 CUMULUS_GRADIENT = vec4(0.01f, 0.0625f, 0.78f, 1.0f);
	float stratus = 1.0f - clamp(cloudType * 2.0f, 0.0, 1.0);
	float stratocumulus = 1.0f - abs(cloudType - 0.5f) * 2.0f;
	float cumulus = clamp(cloudType - 0.5f, 0.0, 1.0) * 2.0f;
	return STRATUS_GRADIENT * stratus + STRATOCUMULUS_GRADIENT * stratocumulus + CUMULUS_GRADIENT * cumulus;
}

float densityHeightGradient(float heightFrac, float cloudType) {
	vec4 cloudGradient = mixGradients(cloudType);
	return smoothstep(cloudGradient.x, cloudGradient.y, heightFrac) - smoothstep(cloudGradient.z, cloudGradient.w, heightFrac);
}

float intersectSphere(vec3 pos, vec3 dir,float r) {
    float a = dot(dir, dir);
    float b = 2.0 * dot(dir, pos);
    float c = dot(pos, pos) - (r * r);
	float d = sqrt((b*b) - 4.0*a*c);
	float p = -b - d;
	float p2 = -b + d;
    return max(p, p2) / (2.0 * a);
}

// Returns density at a given point
// Heavily based on method from Schneider
float density(vec3 pip, vec3 weather, float mip) {
	float time = params.time;
	vec3 p = pip;
	float height_fraction = GetHeightFractionForPoint(length(p));
	p.xz += time * 20.0 * normalize(params.wind_direction) * params.wind_speed * 0.6;
	vec4 n = textureLod(large_scale_noise, p.xyz*0.00008, mip-2.0);
	float fbm = n.g*0.625+n.b*0.25+n.a*0.125;
	float g = densityHeightGradient(height_fraction, weather.r);
	float base_cloud = remap(n.r, -(1.0-fbm), 1.0, 0.0, 1.0);
	float weather_coverage = params.cloud_coverage * weather.b;
	base_cloud = remap(base_cloud*g, 1.0 - (weather_coverage), 1.0, 0.0, 1.0);
	base_cloud *= weather_coverage;
	p.xz -= time * normalize(params.wind_direction) * 40.;
	p.y -= time * 40.;
	vec3 hn = textureLod(small_scale_noise, p*0.001, mip).rgb;
	float hfbm = hn.r*0.625+hn.g*0.25+hn.b*0.125;
	hfbm = mix(hfbm, 1.0-hfbm, clamp(height_fraction*4.0, 0.0, 1.0));
	base_cloud = remap(base_cloud, hfbm*0.4 * height_fraction, 1.0, 0.0, 1.0);
	return pow(clamp(base_cloud, 0.0, 1.0), (1.0 - height_fraction) * 0.8 + 0.5);
}

vec4 march(vec3 pos,  vec3 end, vec3 dir, int depth) {
	const vec3 RANDOM_VECTORS[6] = {vec3( 0.38051305f,  0.92453449f, -0.02111345f),vec3(-0.50625799f, -0.03590792f, -0.86163418f),vec3(-0.32509218f, -0.94557439f,  0.01428793f),vec3( 0.09026238f, -0.27376545f,  0.95755165f),vec3( 0.28128598f,  0.42443639f, -0.86065785f),vec3(-0.16852403f,  0.14748697f,  0.97460106f)};
	float T = 1.0;
	float alpha = 0.0;
	float ss = length(dir);
	dir = normalize(dir);
	vec3 p = pos + dir * hash(pos * 10.0) * ss;
	const float t_dist = sky_t_radius - sky_b_radius;
	float lss = (t_dist / 36.0);
	vec3 ldir = normalize(params.LIGHT_DIRECTION);
	vec3 L = vec3(0.0);
	int count = 0;
	float t = 1.0;
	float costheta = dot(ldir, dir);
	// Stack multiple phase functions to emulate some backscattering
	float phase = max(max(henyey_greenstein(costheta, 0.6), henyey_greenstein(costheta, (0.4 - 1.4 * ldir.y))), henyey_greenstein(costheta, -0.2));
	// Precalculate sun and ambient colors
	// This should really come from a uniform or texture for performance reasons
	vec3 atmosphere_sun = getValFromTLUT(transmittance_lut, tLUTRes, viewPos, params.LIGHT_DIRECTION)*100.0;
	vec3 atmosphere_ambient = atmosphere(normalize(vec3(1.0, 1.0, 0.0))) * 50.0;
	vec3 atmosphere_ground = atmosphere(normalize(vec3(1.0, -1.0, 0.0))) * 50.0;
	
	const float weather_scale = 0.00006;
	float time = params.time * 0.001 + 0.005 * params.time_offset;
	vec2 weather_pos = time * normalize(params.wind_direction) * params.wind_speed;
	
	for (int i = 0; i < depth; i++) {
		p += dir * ss;
		vec3 weather_sample = texture(weather_noise, p.xz * weather_scale + 0.5 + weather_pos).xyz;
		float height_fraction = GetHeightFractionForPoint(length(p));

		t = density(p, weather_sample, 0.0);
		float dt = exp(-params.density * t * ss);
		T *= dt;
		vec3 lp = p;
		float lt = 1.0;
		float cd = 0.0;

		if (t > 0.0) { //calculate lighting, but only when we are in the cloud
			float lheight_fraction = 0.0;
			for (int j = 0; j < 6; j++) {
				lp +=  (ldir + RANDOM_VECTORS[j] * float(j)) * lss;
				lheight_fraction = GetHeightFractionForPoint(length(lp));
				vec3 lweather = texture(weather_noise, lp.xz * weather_scale + 0.5 + weather_pos).xyz;
				lt = density(lp, lweather, float(j));
				cd += lt;
			}
			
			// Take a single distant sample
			lp = p + ldir * 18.0 * lss;
			lheight_fraction = GetHeightFractionForPoint(length(lp));
			vec3 lweather = texture(weather_noise, lp.xz * weather_scale + 0.5).xyz;
			lt = pow(density(lp, lweather, 5.0), (1.0 - lheight_fraction) * 0.8 + 0.5);
			cd += lt;
			
			// captures the direct lighting from the sun
			float beers = exp(-params.density * cd * lss);
			float beers2 = exp(-params.density * cd * lss * 0.25) * 0.7;
			float beers_total = max(beers, beers2);

			vec3 ambient = mix(atmosphere_ground, vec3(1.0), smoothstep(0.0, 1.0, height_fraction)) * params.density * atmosphere_ambient;
			alpha += (1.0 - dt) * (1.0 - alpha);
			L += (ambient + beers_total * atmosphere_sun * phase * alpha) * T * t;
		}
	}
	return clamp(vec4(L, alpha), 0.0, 1.0);
}

// Take a direction as input and draw the sky.
vec4 sky(vec3 dir) {
	vec4 col = vec4(0.0);

	if (dir.y > 0.0) {
		// Only draw clouds above the horizon.
		vec3 camPos = vec3(0.0, g_radius, 0.0);
		vec3 start = camPos + dir * intersectSphere(camPos, dir, sky_b_radius);
		vec3 end = camPos + dir * intersectSphere(camPos, dir, sky_t_radius);
		float shelldist = (length(end - start));
		// Take fewer steps towards horizon
		float steps = (mix(96.0, 54.0, clamp(dot(dir, vec3(0.0, 1.0, 0.0)), 0.0, 1.0)));

		vec3 raystep = dir * shelldist / steps;
		vec4 volume = march(start, end, raystep, int(steps));
		vec3 background = atmosphere(dir);

		// Draw cloud shape
		col = vec4(background * (1.0 - volume.a) + volume.xyz, 1.0);
		// Blend distant clouds into the sky
		col.xyz = mix(clamp(col.xyz, vec3(0.0), vec3(100.0)), clamp(background, vec3(0.0), vec3(100.0)), smoothstep(0.6, 1.3, 1.0 - dir.y));
	} else {
		col = vec4(atmosphere(dir), 1.0);
	}
	
    return col;
}

vec2 oct_wrap(vec2 v) {
	vec2 signVal;
	signVal.x = v.x >= 0.0 ? 1.0 : -1.0;
	signVal.y = v.y >= 0.0 ? 1.0 : -1.0;
	return (1.0 - abs(v.yx)) * signVal;
}

// Hemisphere octahedral. Maximizes use of square texture.
// Adapted from https://johnwhite3d.blogspot.com/2017/10/signed-octahedron-normal-encoding.html
vec3 oct_to_vec3(vec2 e) {
	vec3 n;
	n.x = (e.x - e.y);
	n.y = (e.x + e.y) - 1.0;
	n.z = 1.0 - abs(n.x) - abs(n.y);
    n.xy = n.z >= 0.0 ? n.xy : oct_wrap( n.xy );

	return normalize(n);
}

void main() {
	// Calculate direction from pixel position.
	ivec2 pos = ivec2(gl_GlobalInvocationID.xy) + ivec2(params.update_position);
	vec2 uv = vec2(pos) / params.texture_size;
	vec3 dir = oct_to_vec3(uv).xzy;
	
	imageStore(current_image, pos, sky(dir));

}


